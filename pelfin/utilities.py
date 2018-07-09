import inspect, os, sys, code, traceback
import numpy as np
import json
import csv,re
import importlib, types

RadiiTypes = ['average_all','max_ca_dist','max_heavy_dist']
INF = float('inf')

def recursive_reload(module):
    """Recursively reload modules."""

    reloaded = set()
    def _reload(module):
      try:
        importlib.reload(module)
      except Exception as e:
        print('Warning: {}'.format(e))
      reloaded.add(module)
      for attribute_name in dir(module):
          attribute = getattr(module, attribute_name)
          if type(attribute) is types.ModuleType and attribute not in reloaded:
              _reload(attribute)
    _reload(module)

def gen_pymol_txm(rot, tran):
  """
  Converts BioPython-style rotation and translation into pymol's
  transformation matrix string.

  Args:
  - rot - Bio.PDB.Superimposer().rotran[0]
  - tran - Bio.PDB.Superimposer().rotran[1]

  Returns:
  - _ - string of pymol's transformation matrix.
  """
  rotTp = np.transpose(rot)
  rotTpTran = np.append(rotTp, np.transpose([tran]), axis=1)
  pymolRotMat = np.append(rotTpTran, [[0,0,0,1]], axis=0)
  return '[' + ', '.join(map(str, pymolRotMat.ravel())) + ']'

def int_ceil(f):
  """Ceil a float then turn it into an int."""
  return int(np.ceil(f))

def int_floor(f):
  """Floor a float then turn it into an int."""
  return int(np.floor(f))

def upsample(spec, pts):
  """Upsamples points to be the same number of points in specification. This
  is code translated from Elfin core's C++ code."""
  N = len(spec)

  more_points, fewer_points = (np.copy(spec), np.copy(pts))

  # Compute longer shape total length
  mp_total_length = 0.0
  for i in range(1, N):
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
    while (mp_proportion <= fp_proportion and mpi < N):
      mp_segment = \
        np.linalg.norm(more_points[mpi] - more_points[mpi - 1]) \
        / mp_total_length

      if (mp_proportion + mp_segment) > fp_proportion:
        break
      mp_proportion += mp_segment

      s = (mp_proportion - basefp_proportion) / fp_segment
      upsampled = np.append(upsampled, [base_fp_point + (vec * s)], axis=0)

      mpi += 1

  # Sometimes the last node is automatically added
  if len(upsampled) < N:
    upsampled = np.append(upsampled, [fewer_points[-1]], axis=0)

  return upsampled

def float_approximates(a, b, error=1e-6):
  """Returns whether float a is approximately b within error tolerance"""
  return abs(a-b) < error

def min_dist_from_line(**kwargs):
  """
  Computes the minimum dist from a line to a point.

  Args:
  - point - 3x1 numpy array of the point in quesiton.
  - line_points - Nx(3x1 numpy array) list of the points that define the line
    in quesiton.
  - allow_perp - whether or not the distance can be calculated as a
    perpendicular ray cutting through the lines between points in line_points.
    If false, the distance can only be point-to-point.

  Returns:
  - min_dist - the minimum distance.
  """
  points = kwargs.pop('points')
  line_points = kwargs.pop('line_points')
  allow_perp = kwargs.pop('allow_perp', True)

  min_dist = INF
  for i in range(1, len(line_points)):
    lineSeg = (line_points[i-1], line_points[i])

    # First determine whether point is outside line segment regime
    v = lineSeg[1] - lineSeg[0]
    w = point - lineSeg[0]

    if allow_perp:
      c1 = np.dot(w, v)
      if c1 <= 0: # before lineSeg[0]
        dist = np.linalg.norm(w)
      else:
        c2 = np.dot(v, v)
        if c2 <= c1: # after lineSeg[1]
          dist = np.linalg.norm(point - lineSeg[1])
        else:
          # If not outside, then calculate perpendicular distance
          b = c1 / c2
          pol = lineSeg[0] + b*v
          dist = np.linalg.norm(point - pol)
    else:
      dist = min(np.linalg.norm(point - lineSeg[0]), np.linalg.norm(point - lineSeg[1]))

    if dist < min_dist:
      min_dist = dist

  return min_dist

def check_collision(**kwargs):
  """
  Tests whether a to-be-added node is too close to any node in partially or
  completely formed shape.

  Args:
  - xDB - a dict containing the xDB data. Should have originated from
    read_json().
  - collision_measure - one of RadiiTypes
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

  newCOM = xdb['doubles_data'][nodes[-1]][new_node]['com_b']

  # previous node PAIR (not just single node!) is inherently non-colliding
  for i in range(0, len(nodes) - 2):
    com_dist = np.linalg.norm(shape[i] - newCOM)
    collision_dist = xdb['singles_data'][new_node]['radii'][collision_measure] + \
                        xdb['singles_data'][nodes[i]]['radii'][collision_measure]

    if com_dist < collision_dist:
      return True

  return False

def com_dist_info(xdb):
  """
  Computes centre-of-mass distance information.

  Args:
  - xdb - a dict containing the xdb data. Should have originated from
    read_json().

  Returns:
  - (_,_,_) - tuple containing average, min and max values for centre-of-mass
    distances.
  """
  pd = xdb['doubles_data']
  dists = []
  for s1 in pd.keys():
    for s2 in pd[s1].keys():
      dists.append(np.linalg.norm(pd[s1][s2]['com_b']))

  return np.average(dists), min(dists), max(dists)

def read_csv_points(csvFile):
  """A wrapper of read_csv() but returns as list of numpy array points."""
  pts = []
  
  with open(csvFile, 'r') as file:
    pts = np.asarray([[float(n) for n in re.split(', *| *', l.strip())] for l in file.read().split('\n') if len(l) > 0])
  
  return pts

def read_csv(read_path, delim=','):
  """
  Reads a generic CSV file.

  Args:
  - read_path - string path to read from.
  - delim - delimiter to use for the CSV format.

  Returns:
  - rows - list of rows where each row is a string list of cell values.
  """
  rows = []
  with open(read_path) as csvfile:
    sreader = csv.reader(csvfile, delimiter=delim)
    for r in sreader:
      rows.append([c.strip() for c in r])

  return rows

def save_points_as_csv(**kwargs):
  """
  Saves a list of points into a CSV file.

  Args:
  - points - Nx(3x1 numpy array) list to be saved.
  - save_path - string path to save to.
  - delim - delimiter to use for the CSV format.
  """
  points = kwargs.pop('points')
  save_path = kwargs.pop('save_path')
  delim = kwargs.pop('delim',' ')

  with open(save_path, 'wb') as file:
    wt = csv.writer(file, delimiter=delim)
    for row in points:
      wt.writerow(row)

def read_json(read_path):
  """
  Reads a JSON file.

  Args:
  - read_path - JSON string file path.

  Returns:
  - _ - dict containing data from the JSON file.
  """
  with open(read_path, 'r') as file:
    return json.load(file)

def make_dir(dir):
  """
  Create directory if does not exist.
  """
  if not os.path.exists(dir):
    os.makedirs(dir)

def pause_code(frame=None):
  """
  Pause execution and drop into interactive mode for debugging. This is
  intended to be manually inserted into area of code where debugging is
  needed.

  Args: 
  - frame - specify frame in which the globals and locals are to be debugged.
  """
  print('\n------------------pause_code()------------------')
  if frame is None:
    # Use current frame (one above the exception wrapper)
    frame = inspect.currentframe().f_back
  
  ns = dict(frame.f_globals)
  ns.update(frame.f_locals)
  code.interact(local=ns)

def safe_exec(func, *args):
  """
  Execute func and drops into interactive mode for debugging if an exception
  is raised.

  Args:
  - func - the function handle to be called.
  - *args - args to be expanded for func.
  """
  try:
    func(*args)
  except Exception as e:
    print('\n------------------safe_exec() caught exception------------------')

    # Find last (failed) inner frame
    type, value, tb = sys.exc_info()
    last_frame = lambda tb=tb: last_frame(tb.tb_next) if tb.tb_next else tb
    frame = last_frame().tb_frame
    traceback.print_exc()
    pause_code(frame)

def main():
  raise RuntimeError('This module should not be executed like a script')

if __name__ =='__main__': 
  main()