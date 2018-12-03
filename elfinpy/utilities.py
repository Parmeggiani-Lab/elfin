"""
Elfin data processing utilities module
"""
import inspect
import os
import sys
import code
import traceback as traceback_module
import json
import csv
import re

import numpy as np

RADII_TYPES = ['average_all', 'max_ca_dist', 'max_heavy_dist']
INF = float('inf')

def to_4x4_tx(rot, tran):
    # pause_code()
    tx_3x4 = np.concatenate((rot, np.reshape(tran, (3,1))), axis=1)
    tx_4x4 = np.concatenate((tx_3x4, [[0, 0, 0, 1]]), axis=0)
    return tx_4x4

def dict_diff(A, B):
    if type(A) == list:
        return not all(diff(a, b) for a, b in zip(A, B))
    elif type(A) == dict:
        return not all(diff(A[x], B[x]) for x in A if x in B)
    else:
        return A != B

# https://stackoverflow.com/questions/1036409/recursively-convert-python-object-graph-to-dictionary
def to_dict(obj, classkey=None):
    if isinstance(obj, dict):
        data = {}
        for (k, v) in obj.items():
            data[k] = to_dict(v, classkey)
        return data
    elif hasattr(obj, "_ast"):
        return to_dict(obj._ast())
    elif hasattr(obj, "__iter__") and not isinstance(obj, str):
        return [to_dict(v, classkey) for v in obj]
    elif hasattr(obj, "__dict__"):
        data = dict([(key, to_dict(value, classkey)) 
        for key, value in obj.__dict__.items() 
        if not callable(value) and not key.startswith('_')])
        if classkey is not None and hasattr(obj, "__class__"):
            data[classkey] = obj.__class__.__name__
        return data
    else:
        return obj


def get_rotation(angle_x=0, angle_y=0, angle_z=0):
    """https://en.wikipedia.org/wiki/Rotation_matrix
    """
    radian_x = np.radians(angle_x)
    radian_y = np.radians(angle_y)
    radian_z = np.radians(angle_z)
    rot_x = np.array([
        [1, 0, 0],
        [0, np.cos(radian_x), -np.sin(radian_x)],
        [0, np.sin(radian_x), np.cos(radian_x)]
        ])
    rot_y = np.array([
        [np.cos(radian_y), 0, np.sin(radian_y)],
        [0, 1, 0],
        [-np.sin(radian_y), 0, np.cos(radian_y)]
        ])
    rot_z = np.array([
        [np.cos(radian_z), -np.sin(radian_z), 0],
        [np.sin(radian_z), np.cos(radian_z), 0],
        [0, 0, 1]
        ])

    return np.matmul(a=np.matmul(a=rot_x, b=rot_y), b=rot_z)

def gen_pymol_txm(rot, tran):
    """Converts BioPython-style rotation and translation into pymol's
    transformation matrix string.

    Args:
    - rot - Bio.PDB.Superimposer().rotran[0]
    - tran - Bio.PDB.Superimposer().rotran[1]

    Returns:
    - _ - string of pymol's transformation matrix.
    """
    rot_tp = np.transpose(rot)
    rot_tp_tran = np.append(rot_tp, np.transpose([tran]), axis=1)
    pymol_rot_mat = np.append(rot_tp_tran, [[0, 0, 0, 1]], axis=0)
    return '[' + ', '.join(map(str, pymol_rot_mat.ravel())) + ']'

def int_ceil(float_num):
    """Ceil a float then turn it into an int."""
    return int(np.ceil(float_num))

def int_floor(float_num):
    """Floor a float then turn it into an int."""
    return int(np.floor(float_num))

def upsample(spec, pts):
    """Upsamples points to be the same number of points in specification. This
    is code translated from Elfin core's C++ code.
    """
    n_spec_points = len(spec)

    more_points, fewer_points = (np.copy(spec), np.copy(pts))

    # Compute longer shape total length
    mp_total_length = 0.0
    for i in range(1, n_spec_points):
        mp_total_length += np.linalg.norm(more_points[i] - more_points[i - 1])

    if mp_total_length == INF:
        raise ValueError('Something fishy... mp_total_length is inf!')

    fp_total_length = 0.0
    for i in range(1, len(fewer_points)):
        fp_total_length += np.linalg.norm(fewer_points[i] - fewer_points[i - 1])

    if mp_total_length == INF:
        raise ValueError('Something fishy... fp_total_length is inf!')

    # Upsample fewer_points
    upsampled = np.zeros([0, 3])

    # First and last points are the same
    upsampled = np.append(upsampled, [fewer_points[0]], axis=0)

    mp_proportion = 0.0
    fp_proportion = 0.0
    mpi = 1
    for i in range(1, len(fewer_points)):
        base_fp_point = fewer_points[i - 1]
        next_fp_point = fewer_points[i]
        basefp_proportion = fp_proportion
        fp_segment = np.linalg.norm(next_fp_point - base_fp_point) / fp_total_length
        vec = next_fp_point - base_fp_point

        fp_proportion += fp_segment
        while mp_proportion <= fp_proportion and mpi < n_spec_points:
            mp_segment = \
                np.linalg.norm(more_points[mpi] - more_points[mpi - 1]) \
                / mp_total_length

            if (mp_proportion + mp_segment) > fp_proportion:
                break
            mp_proportion += mp_segment

            scale = (mp_proportion - basefp_proportion) / fp_segment
            upsampled = np.append(upsampled, [base_fp_point + (vec * scale)], axis=0)

            mpi += 1

    # Sometimes the last node is automatically added
    if len(upsampled) < n_spec_points:
        upsampled = np.append(upsampled, [fewer_points[-1]], axis=0)

    return upsampled

def float_approximates(float_a, float_b, error=1e-6):
    """Returns whether float a is approximately b within error tolerance"""
    return abs(float_a-float_b) < error

def check_collision(**kwargs):
    """Tests whether a to-be-added node is too close to any node in partially or
    completely formed shape.

    Args:
    - xDB - a dict containing the xDB data. Should have originated from
        read_json().
    - collision_measure - one of RADII_TYPES
    - nodes - string list of module names
    - new_node - string name the node to be tested
    - shape - Nx(3x1 numpy array) list of node centre-of-masses

    Returns:
    - bool - whether or not the new node, when added to the shape, causes
        collision.
    """
    xdb = kwargs.pop('xdb')
    collision_measure = kwargs.pop('collision_measure')
    nodes = kwargs.pop('nodes')
    new_node = kwargs.pop('new_node')
    shape = kwargs.pop('shape')

    new_com = xdb['double_data'][nodes[-1]][new_node]['com_b']

    # previous node PAIR (not just single node!) is inherently non-colliding
    for i in range(0, len(nodes) - 2):
        com_dist = np.linalg.norm(shape[i] - new_com)
        collision_dist = \
            xdb['single_data'] \
            [new_node]['radii'][collision_measure] + \
            xdb['single_data'] \
            [nodes[i]]['radii'][collision_measure]

        if com_dist < collision_dist:
            return True

    return False

def com_dist_info(xdb):
    """Computes centre-of-mass distance information.

    Args:
    - xdb - a dict containing the xdb data. Should have originated from
        read_json().

    Returns:
    - (_, _, _) - tuple containing average, min and max values for centre-of-mass
        distances.
    """
    double_data = xdb['double_data']
    dists = []
    for single_a_name in double_data.keys():
        for single_b_name in double_data[single_a_name].keys():
            dists.append(
                np.linalg.norm(
                    double_data[single_a_name][single_b_name]\
                        ['com_b']
                    )
                )

    return np.average(dists), min(dists), max(dists)

def read_csv_points(csv_file):
    """A wrapper of read_csv() but returns as list of numpy array points."""
    pts = []

    with open(csv_file, 'r') as file:
        pts = np.asarray(
            [[float(n) for n in re.split(', *| *', l.strip())] \
                for l in file.read().split('\n') if len(l) > 0])

    return pts

def read_csv(read_path, delim=','):
    """Reads a generic CSV file.

    Args:
    - read_path - string path to read from.
    - delim - delimiter to use for the CSV format.

    Returns:
    - rows - list of rows where each row is a string list of cell values.
    """
    rows = []
    with open(read_path) as csv_file:
        sreader = csv.reader(csv_file, delimiter=delim)
        for row in sreader:
            rows.append([c.strip() for c in row])

    return rows

def save_points_as_csv(**kwargs):
    """Saves a list of points into a CSV file.

    Args:
    - points - Nx(3x1 numpy array) list to be saved.
    - save_path - string path to save to.
    - delim - delimiter to use for the CSV format.
    """
    points = kwargs.pop('points')
    save_path = kwargs.pop('save_path')
    delim = kwargs.pop('delim', ' ')

    with open(save_path, 'wb') as file:
        writer = csv.writer(file, delimiter=delim)
        for row in points:
            writer.writerow(row)

def read_json(read_path):
    """Reads a JSON file adn returns a dict."""
    with open(read_path, 'r') as file:
        return json.load(file)

def make_dir(directory):
    """Creates directory if does not exist."""
    if not os.path.exists(directory):
        os.makedirs(directory)

def pause_code(frame=None):
    """Pause execution and drop into interactive mode for debugging. This is
    intended to be manually inserted into area of code where debugging is
    needed.

    Args:
    - frame - specify frame in which the globals and locals are to be debugged.
    """
    print('\n------------------pause_code()------------------')

    if frame is None:
        # Use current frame (one above the exception wrapper)
        frame = inspect.currentframe().f_back

    fi = inspect.getframeinfo(frame)
    print('Where: {loc}:{line}'.format(loc=fi.filename, line=fi.lineno))
    print('What: \n{code}'.format(code=fi.code_context[0]))
    
    name_space = dict(frame.f_globals)
    name_space.update(frame.f_locals)
    code.interact(local=name_space)

def safe_exec(func, *args, **kwargs):
    """Execute func and drops into interactive mode for debugging if an exception
    is raised.

    Args:
    - func - the function handle to be called.
    - *args - args to be expanded for func.
    """
    try:
        func(*args, **kwargs)
    except Exception as ex:
        print('\n------------------safe_exec() caught exception------------------')
        print(ex)

        # Find last (failed) inner frame
        _, _, traceback = sys.exc_info()
        last_frame = \
            traceback.tb_next.tb_next.tb_next \
            if traceback and traceback.tb_next and traceback.tb_next.tb_next \
            else traceback
        frame = last_frame.tb_frame
        traceback_module.print_exc()
        pause_code(frame)

def main():
    """main"""
    raise RuntimeError('This module should not be executed as a script')

if __name__ == '__main__':
    main()
