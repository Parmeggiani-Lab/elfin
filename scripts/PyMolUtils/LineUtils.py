###########################################################
#
#  Pymol script copyright Matthew O'Meara and Xavier Ambroggio 2007
#
#  Last updated Nov 29, 2007
#
#  Draw an axis given a point and a direction.  Optionally give color,
#  length and width.
#
#  Usage: from the pymol menu click file->run...  then find this file.
#  Then at the prompt, type
#
#              draw_axis x,y,z, i,k,j
#
#  where (x,y,z) is the point and (i,k,j) is the direction
#
#
#  Also one can run the script as follows
#
#
#              draw_axis x,y,z, i,k,j, length r,g,b, width,
#
#  where (r,g,b) is the color for the line in red, green, colors,
#  width is the thickness of the line and length is the length of the
#  line.
#
#
#  For a fun example of the script run from the command line after the
#  script is loaded
#
#              draw_axis_example
#
#

from pymol.cgo import *     # get constants
from pymol import cmd
from pymol.vfont import plain

import math
import numpy as np
import json

class Counter:
    def __init__(self):
        self.state = 1
counter = Counter()

def draw(x=None, y=None, z=None, i=None, j=None, k=None, r=1.0, g=1.0, b=1.0, width=1.0, label='', fontSize=10):
    if x == None or y == None or z == None or i == None or j == None or k== None :
        print 'Usage: draw_axis x,y,z, i,k,j, r,g,b, width'
        print 'draw a line centered at (x,y,z) with the direction vector (i,j,k)'
        print 'color (r,g,b), and width arguments are optional'
        # print 'For a fun example of the command, run draw_axis_example'
    else :
        x,y,z = float(x), float(y), float(z)
        i,j,k = float(i), float(j), float(k)
        r,g,b = float(r), float(g), float(b)
        width = float(width)

        # obj = [
        #     LINEWIDTH, width,
        #     BEGIN, LINES,

        #     COLOR,  r,  g,  b,
        #     VERTEX, x, y, z,
        #     VERTEX, i, j, k,

        #     END
        # ]

        w = width
        l = 10
        h = w * 4 # cone hight
        d = w * 2.5 # cone base diameter

        # cone tip
        vStart = np.asarray([x,y,z])
        vEnd = np.asarray([i,j,k])
        dv = vEnd - vStart
        ldv = np.linalg.norm(dv)
        cEnd = vStart + (1+(h/ldv)) * dv

        obj = [
            COLOR,     r,g,b,
            SPHERE,    x,y,z, d, 
            CYLINDER,  x,y,z, i,j,k, w, r,g,b, r,g,b,
            SPHERE,    i,j,k, d, 
        ]

        # # add labels to axes object 
        fontThickness = fontSize / 10
        cyl_text(obj,plain,[i+5,j+5,k],label,
            fontThickness,axes=[[fontSize,0,0],[0,fontSize,0],[0,0,fontSize]])
       
        cmd.load_cgo(obj,'axis'+str(counter.state))
        counter.state += 1

def draw_pts(pts, scale=1.0, width=3.0, color=[0,0,0]):
    pts = np.asarray(pts) * scale
    # delta = 1.0 / len(pts)
    r,g,b = color

    for (p1, p2) in zip(pts, np.roll(pts, -1, axis=0))[0:-1]:
        draw(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], width=width,
            r=r,g=g,b=b)
        # r -= delta
        # b += delta

def draw_csv(specFile, scale=1.0, width=2.0, centred=False, shift=None):
    with open(specFile, 'r') as file:
        pts = np.asarray([[float(n) for n in re.split(', *| *', l.strip())] for l in file.read().split('\n') if len(l) > 0])
        if centred:
            pts = pts - pts[-1]
            dists = [np.linalg.norm(p-[0,0,0]) for p in pts]
            if shift is not None:
                pts += np.asarray(shift) * np.mean(dists)

        draw_pts(pts, scale=scale, width=width)

    cmd.reset()
    cmd.set("depth_cue", 0)

def draw_json(specFile, scale=1.0, width=2.0, centred=False, shift=None):
    with open(specFile, 'r') as file:
        pts = np.asarray(json.load(file)['coms'])
        if centred:
            pts = pts - pts[-1]
            dists = [np.linalg.norm(p-[0,0,0]) for p in pts]
            if shift is not None:
                pts += np.asarray(shift) * np.mean(dists)

        draw_pts(pts, scale=scale, width=width)
        
    cmd.reset()
    cmd.set("depth_cue", 0)

def draw_axis(l=350, w = 1):
    draw(-l,0,0, l,0,0, 1,0,0, w, label='X');
    draw(0,-l,0, 0,l,0, 0,1,0, w, label='Y');
    draw(0,0,-l, 0,0,l, 0,0,1, w, label='Z');


    cmd.reset()
    cmd.set("depth_cue", 0)

def noclip():
    cmd.clip('near', 99999999)
    cmd.clip('far', -99999999)

cmd.extend("draw_axis", draw_axis)
cmd.extend("draw_csv", draw_csv)
cmd.extend("draw_json", draw_json)
draw_axis()
noclip()
print 'LineUtils Loaded'

# a simple example
#draw_line(x=18.232,  y=17.150,  z=9.488,
#             i=-.226639,j=0.708772,k=-.668039,
#             r=1,         b=1,         g=1,
#             width=1,    length=1)


# a more complex example

#import random
#def example1(n, f):
#     """draw a gradient field with n segments with the function f(x,y,z)=(i,j,k)"""
#     for i in range(n):
#          scale = 4
#          x,y,z = [random.random()*scale for i in range(3)]
#          i,j,k = f(x,y,z)

#          draw_axis(x,y,z,i,j,k,abs(i),abs(j),abs(k))


#def f(x,y,z):
#     return (2*x,pow(z,2)+x,y-z)

#cmd.extend("draw_axis_example", lambda :example1(1000,f))
