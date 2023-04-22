#!/usr/bin/env python3
# 
# Test mask for PDK development
# 
# @Authors: Michael Nickerson
# @email: nickersonm@ece.ucsb.edu
# 2022(c)
# 

"""
Test mask for PDK development
(c) Michael Nickerson 2022
"""

## Includes and aliases
# from math import modf
from subprocess import call

import nazca
import marks_NanoFab as marks
import pdk_Fab6 as PDK
cellShift = marks.fCommon.cellShift
cellWidth = marks.fCommon.cellWidth
cellHeight = marks.fCommon.cellHeight


## Definitions
# nazca settings
nazca.cfg.use_hull = True  # Enable convex hull calculations
nazca.font.default_font('nazca')

# Mask-specific definitions
dieName = 'PDK_Test_v1'
dieSize = [[0,0], [2e3, 5e3]]
wgSpace = 25

## Assemble die with assorted components
with nazca.Cell(instantiate=True, name=dieName) as die:
    ## Start adding elements
    y = 500
    
    # Marks
    C = PDK.markVernier
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.markResolution
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.markDEKTAK
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.markTLM
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.markCorner
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    # Interconnects
    y += wgSpace*2
    
    C = PDK.icPassive.strt(500)
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.icBend.strt(500)
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.icIso.strt(500)
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.icShallow.strt(500)
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.icCleave.strt(500)
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.icModulator.strt(500)
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.icActive.strt(500)
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.icActiveDeep.strt(500)
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.icMetalP.strt(500)
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.icMetalTop.strt(500)
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    # Segments
    y += wgSpace
    
    C = PDK.trOpenP()
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.trOpenP(length=20)
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.trPassive2Mod()
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.trShallow2Active
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.segmentIsolation(xs1='wgBend', xs2='wgPassive')
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.segmentIsolation()
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    PDK.icShallow.strt(100).put(0, y)
    C = PDK.segmentSDT().put()
    PDK.icPassive.strt(100).put()
    y += wgSpace*2 + cellHeight(C)
    
    C = PDK.segmentPad()
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    # Components
    y += wgSpace
    
    C = PDK.componentSSC()
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.componentFilter()
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.componentMMI(N=2)
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.componentMMI22()
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.componentMMI(N=3)
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.componentMMI(N=4)
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    PDK.icPassive.strt(50).put(0, y)
    C = PDK.componentUbend().put()
    PDK.icPassive.strt(50).put()
    y += wgSpace + cellHeight(C)
    
    PDK.icPassive.strt(50).put(0, y)
    PDK.componentEcorner().put()
    C = PDK.icPassive.strt(50).put()
    y += wgSpace + cellHeight(C)
    
    PDK.icPassive.strt(50).put(0, y)
    C = PDK.componentCC(xs='wgPassive').put()
    y += wgSpace + cellHeight(C)
    
    C = PDK.componentFan()
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    # Devices
    y += wgSpace
    
    C = PDK.deviceModulator(ploc=0.5)
    C = C.put(0, y + cellHeight(C)/2)
    PDK.segmentInlinePad().put(C.pin['p1'])
    y += wgSpace + cellHeight(C)
    
    C = PDK.deviceModulator()
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.deviceMMI22()
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.deviceMZM()
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.deviceMZM(dual=False)
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.deviceMZM(devicefunc=PDK.deviceSOA)
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.deviceMZM(devicefunc=PDK.deviceSOA, dual=False)
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.devicePMArray(N=3, k=2)
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.devicePMArray(N=2, k=3, soa=500)
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.deviceSOA()
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.deviceSOA(ploc=0)
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.deviceDFB(type='GC', xs2='wgBend')
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.deviceDFB(type='GC', xs2='wgActive')
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.deviceDFB(type='LC')
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.deviceDBR(type='VC')
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.deviceDBR(type='LC')
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.deviceCCDBR(type='VC')
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.deviceCCDBR(type='LC')
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)
    
    C = PDK.deviceCCDBR(xsCC='wgPassive')
    C.put(0, y + cellHeight(C)/2)
    y += wgSpace + cellHeight(C)


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
