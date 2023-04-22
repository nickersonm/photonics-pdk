#!/usr/bin/env python3
# 
# Waveguide test mask for process development
# 
# @Authors: Michael Nickerson
# @email: nickersonm@ece.ucsb.edu
# 2023(c)
# 

"""
Waveguide test mask for process development
    v3 adds CAIBE test mask
    v4 improves circle and star and adds labels
    v5 updates to pdk_Fab4
    v6 shrinks die and removes some items, removes CAIBE test mask, updates resolution block
    v7 updates to pdk_Fab6 and increases resolution
(c) Michael Nickerson 2023
"""

## Includes and aliases
# from math import modf
from subprocess import call

from numpy.core.function_base import linspace

import nazca
import marks_NanoFab as marks
import pdk_Fab6 as PDK
cellShift = marks.fCommon.cellShift
cellWidth = marks.fCommon.cellWidth
cellHeight = marks.fCommon.cellHeight


## Definitions
# nazca settings
# nazca.cfg.use_hull = True  # Enable convex hull calculations
nazca.font.default_font('nazca')

# Mask-specific definitions
dieName = 'WG_Test_v7'
dieSize = [[0,0], [1.75e3, 2.5e3]]

dieLabel = marks.utility.label(text=dieName, height=50, layer=1, date=True)

# Waveguide test definitions
minWidth = 0
nomWidth = 2
maxWidth = 4
ringSize = 750
taperR0 = 0
taperLen = ringSize/2 - taperR0
wgLength = 825


## Simple parts
# Tapered waveguide
with nazca.Cell(name='Taper: Linear', instantiate=True) as taperLinear:
    PDK.icPassive.taper(length=taperLen, width1=minWidth, width2=maxWidth).put(0)
    # Put width mark locations
    for w in [0.5, 1, 1.5, 2, 4]:
        x = taperLen*(w - minWidth)/(maxWidth - minWidth)
        PDK.icPassive.strt(length=5, width=1, arrow=False).put(x, 3, 90)
        PDK.icPassive.strt(length=5, width=1, arrow=False).put(x, -3, -90)

# Ring
with nazca.Cell(name='Ring', instantiate=True) as ring:
    marks.fCommon.layerPolygon(1, nazca.geom.ring(radius=ringSize/2, width=nomWidth*2, N=256)).put(0)


## Multi-component parts
# MMI to array with fan, various kinds
with nazca.Cell(name='Waveguides: Fan') as groupFan:
    y=-25
    PDK.icBend.strt(50, arrow=False).put(0, y)
    m = PDK.deviceMMI_tree(N=2, k=2, pitch=10, straight=0).put()
    f = PDK.componentFan(N=2**2, pitchin=10, pitchout=3).put()
    for i in range(2**2):
        PDK.straightToOut(wgLength, node=f.pin['b'+str(i)], xs='wgBend', ssc=False).put(f.pin['b'+str(i)])
    y -= 150
    
    PDK.icBend.strt(50, arrow=False).put(0, y)
    m = PDK.deviceMMI_tree(N=2, k=3, pitch=10, straight=0).put()
    f = PDK.componentFan(N=2**3, pitchin=10, pitchout=4).put()
    for i in range(2**3):
        PDK.straightToOut(wgLength, node=f.pin['b'+str(i)], xs='wgBend', ssc=False).put(f.pin['b'+str(i)])
    y -= 150
    
    PDK.icBend.strt(50, arrow=False).put(0, y)
    m = PDK.deviceMMI_tree(N=4, k=1, pitch=10, straight=0).put()
    f = PDK.componentFan(N=4, pitchin=10, pitchout=4).put()
    for i in range(4):
        PDK.straightToOut(wgLength, node=f.pin['b'+str(i)], xs='wgBend', ssc=False).put(f.pin['b'+str(i)])
    y -= 150
    
    PDK.icBend.strt(50, arrow=False).put(0, y)
    m = PDK.deviceMMI_tree(N=3, k=2, pitch=10, straight=0).put()
    f = PDK.componentFan(N=3**2, pitchin=10, pitchout=3).put()
    for i in range(3**2):
        PDK.straightToOut(wgLength, node=f.pin['b'+str(i)], xs='wgBend', ssc=False).put(f.pin['b'+str(i)])
groupFan = groupFan.remove_layer(range(2,15))

# Assorted s-bends
with nazca.Cell(name='Waveguides: S-bend') as groupSbend:
    y = 0
    for br in range(10, 200, 10):
        b = PDK.icPassive.sinebend(distance=PDK.sinebend_minlen(100, br), offset=100, arrow=False)
        b1 = b.put(wgLength/2, y-100, flop=True)
        PDK.icPassive.strt(b1.pinout.x, arrow=False).put(b1.pinout)
        b1 = b.put(b1.pinin)
        
        PDK.icPassive.strt(wgLength - b1.pinout.x, arrow=False).put(b1.pinout)
        y -= 25
groupSbend = groupSbend.remove_layer(range(2,15))

# Assorted Curves
with nazca.Cell(name='Waveguides: Curves') as groupCurves:
    y = 0
    for br in range(300, 50, -25):
        PDK.icPassive.strt(wgLength - br - 75, arrow=False).put(0, y)
        PDK.icPassive.euler_bend(angle=90, radius=br, instantiate=False, arrow=False).put()
        PDK.icPassive.strt(dieSize[1][0] - 100 - nazca.cp.y(), arrow=False).put()
        y += 50
groupCurves = groupCurves.remove_layer(range(2,15))

## Assemble die
with nazca.Cell(instantiate=True, name=dieName) as die:
    ## Utility marks
    # Die outline
    marks.utility.die(dieSize).put(0,0)
    marks.utility.die(dieSize, layer=1005, grow=300).put(0,0)
    
    # Corner marks
    marks.fCommon.putCorners(PDK.markCorner.remove_layer([15]), 
                             inset=0, flip=True, flop=True, diesize=dieSize)
    
    # Die label on top and bottom
    dieLabel.put(dieSize[1][0]/2, 50)
    dieLabel.put(dieSize[1][0]/2, dieSize[1][1] - 50 - cellHeight(dieLabel))
    
    # Alignment marks in explicit locations
    for x in [dieSize[0][0]+100, dieSize[1][0]-100]:
        for y in [dieSize[0][1]+100, dieSize[1][1]-100]:
            marks.mla150.AlignHighMag().put(x,y)
    
    # Resolution test and DEKTAK in all 4 corners
    marks.fCommon.putCorners(
        marks.utility.ResolutionBlockSet(xs_pattern=1, res=[0.6, 0.8, 0.9, 1, 1.1, 1.3]), 
        inset=[250, 25], shift=True, diesize=dieSize)
    marks.fCommon.putCorners(
        marks.utility.DEKTAK_box(xs_pad=1), 
        inset=[400, 25], shift=True, diesize=dieSize)
    
    
    ## Elements
    # Concentric rings
    for s in linspace(0.1, 1, 10):
        ring.put(dieSize[1][0]/4, dieSize[1][1] - ringSize/2 - 250, scale=s)
        marks.utility.label('%.2fµm' % (s*nomWidth*2), height=20, date=False).\
            put(dieSize[1][0]/4, dieSize[1][1] - ringSize/2 - 250 + s*ringSize/2 + 5)
    
    # Radiating linear tapers
    for i in range(16):
        nazca.cp.goto(dieSize[1][0]*3/4, dieSize[1][1] - ringSize/2 - 250, i*45/2)
        nazca.cp.skip(taperR0)
        taperLinear.put()
    # Mark width locations
    for w in [0.5, 1, 1.5, 2, 4]:
        r = taperR0 + taperLen*(w - minWidth)/(maxWidth - minWidth)
        cellShift(marks.utility.label('%.2fµm' % w, height=10, date=False), ['left', 'lower']).\
            put(dieSize[1][0]*3/4 + r + 5, dieSize[1][1] - ringSize/2 - 250 + 5)
    
    # S-bends
    groupSbend.put(50, 200, 90)
    
    # MMI Fan
    groupFan.put(650, 200, 90)
    
    # Curves
    groupCurves.put(dieSize[1][0]-50, 200, 90)
    x += cellHeight(groupCurves)
die = cellShift(die, 'center')  # Center die for easier testing


## Generate GDS
from datetime import datetime
dieFile = dieName + '_' + datetime.now().strftime('%Y%m%d') + '.gds'
nazca.export_gds(topcells=die, filename=dieFile, clear=False)

# Postprocess with KLayout
print('Postprocessing with KLayout...')
from subprocess import call
from os.path import abspath
postDRC = [abspath('..\\pdk_Fab6\\Fab6_postprocess.lydrc'),
           abspath('..\\pdk_Fab6\\Fab6_postprocess_merge.lydrc')]
for drc in postDRC:
    call(['klayout_app.exe', '-b', '-r', drc, '-rd', 'input='+abspath(dieFile)])
