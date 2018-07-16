#!/usr/bin/env python3

#
# A PyMol extension script for drawing lines.
#

def main():
    """main"""
    raise RuntimeError('This module should not be executed as a script')

if __name__ =='__main__': 
    main()

in_pymol = False
try:
    import pymol
    in_pymol = True
except ImportError as ie:
    main()

if in_pymol:
    import pymol.cgo # constants
    from pymol import cmd
    from pymol.vfont import plain

    import numpy as np
    import json

    class LineCounter:
        """A counter to keep track of how many lines are drawn and to be used in naming
        the next line.
        """
        def __init__(self):
            self.count = 1
        def add(self):
            self.count += 1
        def get(self):
            return self.count
    counter = LineCounter()

    @cmd.extend
    def draw_line(starting_point=None, line_vector=None, color=(1,1,1), width=1.0, label='', font_size=12):
        """Draws a line.

        Args: 
        - starting_point - 3-value list or tuple (x,y,z)
        - line_vector - 3-value list or tuple (i,j,k)
        - color - 3-value list or tuple (r,g,b)
        - width - float
        - label - stirng
        - font_size - float
        """
        if starting_point is None or len(starting_point) != 3 \
            or line_vector is None or len(line_vector) != 3:
            print(draw_line.__doc__)
        else :
            x,y,z = starting_point
            i,j,k = line_vector
            r,g,b = color

            d = width * 2.5 # cone base diameter

            obj = [
                pymol.cgo.COLOR,     r,g,b,
                pymol.cgo.SPHERE,    x,y,z, d, 
                pymol.cgo.CYLINDER,  x,y,z, i,j,k, width, r,g,b, r,g,b,
                pymol.cgo.SPHERE,    i,j,k, d, 
            ]

            # # add labels to axes object 
            font_thickness = font_size / 10
            pymol.cgo.cyl_text(obj,plain,[i+5,j+5,k],label,
                font_thickness,axes=[[font_size,0,0],[0,font_size,0],[0,0,font_size]])
             
            cmd.load_cgo(obj,'line'+str(counter.get()))
            counter.add()

    @cmd.extend
    def draw_points(points=[], scale=1.0, width=3.0, color=(0,0,0)):
        """Draws points and joins them up using colored lines.

        Args:
        - points - list of 3-value lists
        - scale - float
        - width - float
        - color - 3-value list or tuple
        """
        points = np.asarray(points) * scale
        r,g,b = color

        for (p1, p2) in zip(points, np.roll(points, -1, axis=0))[0:-1]:
            draw_line(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], width=width,
                r=r,g=g,b=b)

    @cmd.extend
    def draw_csv(spec_file=None, scale=1.0, width=2.0, centered=False, shift=None):
        """Draws points specified by a csv.

        Args:
        - spec_file - string path
        - scale - float
        - width - float
        - centered - boolean
        - shift - 3-value list or tuple
        """
        if spec_file is None:
            print(draw_csv.__doc__)
        else:
            with open(spec_file, 'r') as file:
                pts = np.asarray([[float(n) for n in re.split(', *| *', l.strip())] for l in file.read().split('\n') if len(l) > 0])
                if centered:
                    pts = pts - pts[-1]
                    dists = [np.linalg.norm(p-[0,0,0]) for p in pts]
                    if shift is not None:
                        pts += np.asarray(shift) * np.mean(dists)

                draw_points(pts, scale=scale, width=width)

            cmd.reset()
            cmd.set("depth_cue", 0)

    @cmd.extend
    def draw_json(spec_file, scale=1.0, width=2.0, centered=False, shift=None):
        """Draws points specified by a json.

        Args:
        - spec_file - string path
        - scale - float
        - width - float
        - centered - boolean
        - shift - 3-value list or tuple
        """
        if spec_file is None:
            print(draw_csv.__doc__)
        else:
            with open(spec_file, 'r') as file:
                pts = np.asarray(json.load(file)['coms'])
                if centered:
                    pts = pts - pts[-1]
                    dists = [np.linalg.norm(p-[0,0,0]) for p in pts]
                    if shift is not None:
                        pts += np.asarray(shift) * np.mean(dists)

                draw_points(pts, scale=scale, width=width)
                
            cmd.reset()
            cmd.set("depth_cue", 0)

    @cmd.extend
    def draw_axes(length=500, width=2, font_size=20):
        """Draws the XYZ axes."""
        draw_line(
            starting_point=(-length,0,0), 
            line_vector=(length,0,0), 
            color=(1,0,0), 
            width=width, 
            label='X', 
            font_size=font_size);
        draw_line(
            starting_point=(0,-length,0), 
            line_vector=(0,length,0), 
            color=(0,1,0), 
            width=width, 
            label='Y', 
            font_size=font_size);
        draw_line(
            starting_point=(0,0,-length), 
            line_vector=(0,0,length), 
            color=(0,0,1), 
            width=width, 
            label='Z', 
            font_size=font_size);

        cmd.reset()
        cmd.set("depth_cue", 0)

    @cmd.extend
    def noclip():
        """Sets clipping to nearly infinity."""
        cmd.clip('near', 99999999)
        cmd.clip('far', -99999999)

    draw_axes()
    noclip()

    print('Line Utils Extension Loaded')