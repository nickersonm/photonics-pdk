#!/usr/bin/env python3
# 
# Fab6 laser, SOA, and MZM test structures
# 
# @Authors: Michael Nickerson
# @email: nickersonm@ece.ucsb.edu
# 2023(c)
# 

"""
Fab6 laser, SOA, and MZM test structures for "AR8.2" GaAs epitaxy
(c) Michael Nickerson 2023
    v2: Simplified test structures and reduced to optimal laser types
    v3: Added ring laser
"""

## Includes and aliases
from datetime import datetime
from math import ceil

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
dieName = 'Active_Test_v3'
dieSize = [[0,0], [3e3, 2.65e3]]
wgSpace = PDK.contactWidth + PDK.metalBuffer # Waveguide spacing
wgOverlap = PDK.cleaveWidth*1.5  # Overlap between subdies, 1/2 space between cleaving lines, and general inset
activeLen = 750 # Typical active length
componentSpace = 100    # Typical longitudinal distance between components
angleOut = 6    # Output facet angle to avoid backreflections

# Derived definitions
dieLabel = marks.utility.label(text=dieName, height=75, layer='MetalTop', date=True, origin=['top', 'left'])



### Die subsections
## Test section for single-sided AR devices, namely lasers
with nazca.Cell(name='Laser Test Devices') as groupTestLaser:
    subLen = dieSize[1][0]  # Full die length
    d, iD = 0, 0
    x, y = 75, wgSpace+5    # Extra y space to compensate for angled outputs
    
    for lg in [0.5, 0.8]:
        # VC-DBRs
        PDK.inputLabel(text='%i.%02i: %.3g mm VC-DBR / %.3g mm active' % (d, iD, lg, 1e-3*activeLen), 
                       xs=None, texty=PDK.metalBuffer*2).put(x, y)
        PDK.icPassive.strt(length=wgOverlap*2).put(x, y)
        PDK.deviceDBR(type='VC', gratinglength=lg*1e3, activelength=activeLen).put()
        PDK.soaToOut(x=subLen, soalen=activeLen/2, angle=angleOut).put()
        iD += 1
        y += 2*wgSpace
        
        # LC-DBRs
        PDK.inputLabel(text='%i.%02i: %.3g mm LC-DBR / %.3g mm active' % (d, iD, lg, 1e-3*activeLen), 
                       xs=None, texty=PDK.metalBuffer*2).put(x, y)
        PDK.icShallow.strt(length=wgOverlap*2).put(x, y)
        PDK.deviceDBR(type='LC', gratinglength=lg*1e3, activelength=activeLen).put()
        PDK.segmentSDT().put()
        PDK.soaToOut(x=subLen, soalen=activeLen/2, angle=angleOut).put()
        iD += 1
        y += 2*wgSpace
        
    for lg in [1.0, 1.5]:
        # GC-DFBs
        PDK.inputLabel(text='%i.%02i: %.3g mm GC-DFB' % (d, iD, lg), 
                       xs=None, texty=PDK.metalBuffer*2).put(x, y)
        PDK.icShallow.strt(length=wgOverlap*2).put(x, y)
        PDK.deviceDFB(length=lg*1e3).put()
        PDK.soaToOut(x=subLen, soalen=activeLen/2, angle=angleOut).put()
        iD += 1
        y += 2*wgSpace
    
    # CC lasers
    lg = 0.30
    PDK.inputLabel(text='%i.%02i: CC / %.3g mm active / %.3g mm VC-DBR' % (d, iD, 1e-3*activeLen, lg), 
                   xs=None, texty=PDK.metalBuffer*2).put(x, y)
    PDK.deviceCCDBR(type='VC', gratinglength=lg*1e3, activelength=activeLen).put(x, y)
    PDK.soaToOut(x=subLen, soalen=activeLen/2, angle=angleOut).put()
    iD += 1
    y += 2*wgSpace
    
    lg = 0.50
    PDK.inputLabel(text='%i.%02i: CC / %.3g mm active / %.3g mm LC-DBR' % (d, iD, 1e-3*activeLen, lg), 
                   xs=None, texty=PDK.metalBuffer*2).put(x, y)
    PDK.deviceCCDBR(type='LC', gratinglength=lg*1e3, activelength=activeLen).put(x, y)
    PDK.segmentSDT().put()
    PDK.soaToOut(x=subLen, soalen=activeLen/2, angle=angleOut).put()
    iD += 1
    y += 2*wgSpace
    
    # Restricted-regrowth CC laser
    lg = 0.30
    PDK.inputLabel(text='%i.%02i: CC / %.3g mm active-restricted / %.3g mm VC-DBR' % (d, iD, 1e-3*activeLen, lg), 
                   xs=None, texty=PDK.metalBuffer*2).put(x, y)
    L = PDK.deviceCCDBR(type='VC', gratinglength=lg*1e3, activelength=activeLen, ploc=0).\
        remove_layer('ProtectRegrowth').put(x, y)
    PDK.soaToOut(x=subLen, soalen=activeLen/2, angle=angleOut).put()
    marks.utility.layerPolygon(
        poly=nazca.geometries.trapezoid(length=20, 
                                        height=4, 
                                        position=8, angle1=60, angle2=90), 
        layers='ProtectRegrowth').put(L.pin['p0'].x + 20, y)
    marks.utility.layerPolygon(
        poly=nazca.geometries.rectangle(length=L.pin['p1'].x-L.pin['p0'].x-1.5*20, 
                                        height=4, position=2), 
        layers='ProtectRegrowth').put(L.pin['p0'].x + 20, y)
    marks.utility.layerPolygon(
        poly=nazca.geometries.trapezoid(length=20, 
                                        height=4, 
                                        position=8, angle1=90, angle2=60), 
        layers='ProtectRegrowth').put(L.pin['p1'].x, y)
    iD += 1
    y += 2*wgSpace
    
    # Etched-facet lasers
    PDK.inputLabel(text='%i.%02i: Etch / %.3g mm active / %.3g mm VC-DBR' % (d, iD, 1e-3*activeLen, lg), 
                   xs=None, texty=PDK.metalBuffer*2).put(x, y)
    PDK.icShallow.strt(length=4*wgOverlap).put(x, y)
    PDK.deviceDBR(type='VC', rfrac=1, gratinglength=lg*1e3, activelength=activeLen).put()
    PDK.soaToOut(x=subLen, soalen=activeLen/2, angle=angleOut).put()
    iD += 1
    y += 2*wgSpace
    
    # Cleaved-facet lasers
    PDK.inputLabel(text='%i.%02i: Cleave / %.3g mm active / %.3g mm VC-DBR' % (d, iD, 1e-3*activeLen, lg), 
                   xs=None, ssc=False, texty=PDK.metalBuffer*2).put(x, y)
    PDK.icShallow.strt(length=4*wgOverlap).put(0, y)
    PDK.deviceDBR(type='VC', rfrac=1, gratinglength=lg*1e3, activelength=activeLen).put(x, y)
    PDK.soaToOut(x=subLen, soalen=activeLen/2, angle=angleOut).put()
    iD += 1
    y += 2*wgSpace
    
    # Deep ridge ring laser
    PDK.inputLabel(text='%i.%02i: %.3g mm active ring laser' % (d, iD, 1e-3*(2*1500)), 
                   xs=None, texty=PDK.metalBuffer*2).put(x, y)
    with nazca.Cell(name="RingLaser", instantiate=True) as devRingLaser:
        # Initial U-bend
        U = PDK.componentUbend(width=3, Lout=4).put(0, 11, flop=True)
        
        # Straight sections with transitions
        PDK.trPassive2ActiveDeep(w1=3, length=5).put(U.pinout)
        PDK.icActiveDeep.strt(length=1500, arrow=False).put()
        R = PDK.trPassive2ActiveDeep(xs1='wgPassive', length=5).flip().put()
        PDK.trPassive2ActiveDeep(w1=3, length=5).put(U.pinin)
        PDK.icActiveDeep.strt(length=1500, arrow=False).put()
        PDK.trPassive2ActiveDeep(xs1='wgBend', length=5).flip().put()
        
        # MMI split
        bendOff = 1.5
        bendLen = ceil(PDK.sinebend_minlen(dy=bendOff, radius=PDK.wgBend.radius))
        M = PDK.componentMMI(N=2, xs='wgBend').put()
        PDK.icPassive.sinebend_dw(length=bendLen, offset=bendOff, 
                                  width1=1.5, width2=2).put(M.pin['b1']).raise_pins(['b0'])
        PDK.icPassive.sinebend_dw(length=bendLen, offset=-bendOff, 
                                  width1=1.5, width2=2).put(M.pin['b0'])
        
        # Loop back to other arm
        PDK.componentUbend(width=2, Lout=2).put()
        PDK.icPassive.sbend_p2p(pin2=R.pinout, arrow=False).put()
        
    devRingLaser.put(140, y)
    PDK.soaToOut(x=subLen, soalen=activeLen/2, angle=angleOut).put()
    iD += 1
    y += 2*wgSpace


## Test section for double-sided AR devices
with nazca.Cell(name='Dual-AR Test Devices') as groupTestAR:
    subLen = dieSize[1][0]  # Full die length
    d, iD = 1, 0
    x, y = 0, wgSpace
    
    # MZMs
    for lg in [0.5, 1, 1.8]:
        PDK.inputLabel(text='%i.%02i: %.3g mm MZM' % (d, iD, lg), textx=componentSpace, 
                       ssc=True, straight=componentSpace).put(x, y)
        PDK.componentFilter().put()
        PDK.icPassive.strt(length=20, arrow=False).put()
        PDK.deviceMZM(length=lg*1e3, dual=False).put()
        PDK.straightToOut(x=subLen).put()
        iD += 1
        y += 2*wgSpace
    
    # Straight SOAs for crude amplitude measurement
    for lg in [0.5, 1, 1.5, 2]:
        PDK.inputLabel(text='%i.%02i: %.3g mm SOA' % (d, iD, lg), textx=componentSpace, 
                       ssc=True, straight=componentSpace).put(x, y)
        PDK.componentFilter().put()
        PDK.icPassive.strt(length=componentSpace, arrow=False).put()
        PDK.deviceSOA(length=lg*1e3).put()
        PDK.straightToOut(x=subLen).put()
        iD += 1
        y += 2*wgSpace
    
    # Restricted-regrowth SOA
    lg = 1.5
    PDK.inputLabel(text='%i.%02i: %.3g mm restricted-MQW SOA' % (d, iD, lg), textx=componentSpace, 
                    ssc=True, straight=componentSpace).put(x, y)
    PDK.componentFilter().put()
    PDK.icPassive.strt(length=componentSpace, arrow=False).put()
    L = PDK.deviceSOA(length=lg*1e3, ploc=0).remove_layer('ProtectRegrowth').put()
    PDK.straightToOut(x=subLen).put()
    marks.utility.layerPolygon(
        poly=nazca.geometries.trapezoid(length=20, 
                                        height=4, 
                                        position=8, angle1=60, angle2=90), 
        layers='ProtectRegrowth').put(L.pin['p0'].x + 20, y)
    marks.utility.layerPolygon(
        poly=nazca.geometries.rectangle(length=L.pin['p1'].x-L.pin['p0'].x-1.5*20, 
                                        height=4, position=2), 
        layers='ProtectRegrowth').put(L.pin['p0'].x + 20, y)
    marks.utility.layerPolygon(
        poly=nazca.geometries.trapezoid(length=20, 
                                        height=4, 
                                        position=8, angle1=90, angle2=60), 
        layers='ProtectRegrowth').put(L.pin['p1'].x, y)
    iD += 1
    y += 2*wgSpace
    
    # MZM-SOAs for SOA phase measurement
    for lg in [1.0, 1.5]:
        PDK.inputLabel(text='%i.%02i: %.3g mm MZM-SOA' % (d, iD, lg), textx=componentSpace, 
                       ssc=True, straight=componentSpace).put(x, y)
        PDK.componentFilter().put()
        PDK.icPassive.strt(length=20, arrow=False).put()
        PDK.deviceMZM(length=lg*1e3, dual=False, devicefunc=PDK.deviceSOA).put()
        PDK.straightToOut(x=subLen).put()
        iD += 1
        y += 2*wgSpace
    
    # MZM + MZ-SOA
    lg = 0.5
    PDK.inputLabel(text='%i.%02i: %.3g mm MZM + %.3g mm SOA' % (d, iD, lg, lg), textx=componentSpace, 
                    ssc=True, straight=componentSpace*3, xs='wgBend').put(x, y)
    PDK.componentFilter(xs='wgBend').put()
    PDK.icBend.strt(length=20, arrow=False).put()
    M = PDK.deviceMMI(N=2, xs='wgBend').put()
    P = PDK.deviceModulator(length=lg*1e3, xs1='wgBend')
    P.put(M.pin['b1'])
    PDK.icBend.strt(length=P.geolen()).put(M.pin['b0'])
    M = PDK.deviceMMI22(xs='wgBend').put()
    P = PDK.deviceSOA(length=lg*1e3, xs1='wgBend')
    P.put(M.pin['b1'])
    PDK.icBend.strt(length=P.geolen()).put(M.pin['b0'])
    PDK.deviceMMI(N=2, xs='wgBend').flip().put(flip=True)
    PDK.straightToOut(x=subLen, xs='wgBend').put()
    iD += 1
    y += 2*wgSpace



### Assemble die
with nazca.Cell(instantiate=True, name=dieName) as die:
    ## Utility marks
    # Die outline for lightfield and darkfield
    marks.utility.die(dieSize).put(0,0)
    marks.utility.die(dieSize, layer=1005, grow=200).put(0,0)
    
    # MLA150 alignment marks in explicit asymmetric locations
    for x, y in [[200, 75], [dieSize[1][0]-125, 75],
                 [100, dieSize[1][1]-85], [dieSize[1][0]-110, dieSize[1][1]-70]]:
        marks.mla150.AlignHighMag(layer=['ProtectRegrowth', 'ProtectRidge']).put(x,y)
        nazca.netlist.Annotation(layer='Annotation', 
                                    text=( "Align: %.0f, %.0f" % (x,y) )).put(x,y)
    
    # Stepper2 DFAS marks in explicit asymmetric locations
    for x, y in [[1800, 50], [1950, dieSize[1][1]-50]]:
        marks.stepper2.LocalAlign(layer=['ProtectRidge']).put(x,y)
        nazca.netlist.Annotation(layer='Annotation', 
                                    text=( "Align: %.0f, %.0f" % (x,y) )).put(x,y)
    
    y = 15
    # Resolution tests and verniers in lower left and right
    marks.fCommon.putCorners(PDK.markResolution, diesize=dieSize, inset=[300, y], quadrants=[3])
    marks.fCommon.putCorners(PDK.markVernier, diesize=dieSize, inset=[225, y], quadrants=[4])
    
    # Layer labels on the top left
    marks.fCommon.putCorners(marks.utility.LayerLabels(layers=[1,15,2]),
                             diesize=dieSize, inset=[250, y], quadrants=[2])
    marks.fCommon.putCorners(marks.utility.LayerLabels(layers=[3,4,5]),
                             diesize=dieSize, inset=[525, y], quadrants=[2])
    
    # Custom TLMs in lower middle and upper right
    # On ridge
    marks.fCommon.putCorners(marks.utility.ConcentricTLM(xs_pad=['MetalTop', 'MetalVia'], 
                                                         grow_pad=[2, 0], 
                                                         xs_background='ProtectRidge',
                                                         grow_background=10), 
                             diesize=dieSize, inset=[1250, -40], quadrants=[3])
    marks.utility.putCorners(marks.utility.CircularTLM(xs_pad=['MetalTop', 'MetalVia'], 
                                                       grow_pad=[1, 0], buffer=40, 
                                                       xs_background='ProtectRidge', 
                                                       grow_background=4, 
                                                       r0=[8], ratio=[1, 1.5, 2.0],
                                                       ratio_outer=1.5), 
                             diesize=dieSize, inset=[250, 9], quadrants=[1])
    marks.utility.putCorners(marks.utility.CircularTLM(xs_pad=['MetalTop', 'MetalVia'], 
                                                       grow_pad=[1, 0], buffer=40, 
                                                       xs_background='ProtectRidge', 
                                                       grow_background=4, 
                                                       r0=[14], ratio=[1, 1.5, 2.0],
                                                       ratio_outer=1.5), 
                             diesize=dieSize, inset=[250+325, 0], quadrants=[1])
    # Off ridge
    marks.fCommon.putCorners(marks.utility.ConcentricTLM(xs_pad=['MetalTop', 'MetalVia'], 
                                                         grow_pad=[2, 0]), 
                             diesize=dieSize, inset=[1500, -35], quadrants=[3])
    marks.utility.putCorners(marks.utility.CircularTLM(xs_pad=['MetalTop', 'MetalVia'], 
                                                       grow_pad=[1, 0], buffer=40, 
                                                       r0=[14], ratio=[1, 1.5, 2.0],
                                                       ratio_outer=1.5), 
                             diesize=dieSize, inset=[750, 0], quadrants=[4])
    
    # Die label in upper middle
    dieLabel.put(800, dieSize[1][1] - y)
    
    # Horizontal cleave lanes for test marking separation
    y = 175
    PDK.icCleave.strt(length=dieSize[1][0] + 2*wgOverlap, arrow=False).put(-wgOverlap, y)
    PDK.icCleave.strt(length=dieSize[1][0] + 2*wgOverlap, arrow=False).put(-wgOverlap, 2475 - wgOverlap/2)
    y += wgSpace
    
    # Die border facet cleave lanes
    PDK.icCleave.strt(length=dieSize[1][1] + wgOverlap*2, arrow=False).put(wgOverlap/2, -wgOverlap, 90)
    PDK.icCleave.strt(length=dieSize[1][1] + wgOverlap*2, arrow=False).put(dieSize[1][0] - wgOverlap/2, -wgOverlap, 90)
    
    
    ## Devices
    # Laser test devices
    groupTestLaser.put(0, y)
    y += 1150
    
    # Optional laser facet break cleave
    PDK.icCleave.strt(length=y + 2*wgOverlap, arrow=False).put(75 + wgOverlap, y + wgOverlap, -90)
    
    # Separate with cleave line
    PDK.icCleave.strt(length=dieSize[1][0] + wgOverlap*2, arrow=False).put(-wgOverlap, y, 0)
    y += wgSpace
    
    # MZM & SOA test devices
    groupTestAR.put(0, y)
    y += cellHeight(groupTestAR)


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
