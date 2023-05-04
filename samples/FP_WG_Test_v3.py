#!/usr/bin/env python3
# 
# FP waveguide test for phase modulation coefficient extraction
# 
# @Authors: Michael Nickerson
# @email: nickersonm@ece.ucsb.edu
# 2023(c)
# 

"""
FP waveguide test for phase modulation coefficient extraction
(c) Michael Nickerson 2023
    v3: remove isolation etches, add SSC to all but SSC=0 test, constant SDT shallow length, added Ubend
"""

## Includes and aliases
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
dieName = 'FP_WG_Test_v3'
dieSize = [[0,0], [None, 2.50e3]]   # Simple FP waveguide tests
dieSize[1][0] = dieSize[1][1]*2 - 200   # Adjust for test markings
modLen = 1000   # Typical modulator length
wgSpace = PDK.contactWidth + PDK.metalBuffer # Waveguide spacing
wgOverlap = PDK.cleaveWidth*1.5  # Overlap between subdies, 1/2 space between cleaving lines, and general inset



### Die subsections
## FP structures to measure loss
with nazca.Cell('FP Tests') as groupFP:
    subLen = dieSize[1][1]  # Subdie length = die height to accomodate rotation
    d, i = 2, 0
    x, y, straight = 0, wgSpace, 100
    
    # SSC
    for b in [True, False]:
        PDK.inputLabel(text='SSC: %i' % (b*2), ssc=b, straight=straight).put(0, y)
        PDK.modToOut(x=subLen, ssc=b, modlen=modLen, isolation=False).put()
        y += wgSpace*2
    
    # Mode filter
    for n in [1, 3, 6]:
        PDK.inputLabel(text='Mode Filter: %i' % n, straight=straight).put(0, y)
        for _ in range(n):
            PDK.componentFilter().put()
            PDK.icPassive.strt(length=wgOverlap).put()
        PDK.modToOut(x=subLen, modlen=modLen, isolation=False).put()
        y += wgSpace*2
    
    # Bend radius, 2 µm WG
    for br in [200, 100, 50, 25]:
        PDK.inputLabel(text='2 µm WG: %.0f µm BR' % br, straight=straight).put(0, y)
        PDK.icPassive.strt(300 - PDK.sinebend_minlen(wgSpace, br), arrow=False).put()
        PDK.icPassive.sinebend(distance=PDK.sinebend_minlen(wgSpace, br), offset=-wgSpace, arrow=False).put()
        PDK.icPassive.sinebend(distance=PDK.sinebend_minlen(wgSpace, br), offset=-wgSpace, arrow=False).put(flip=True)
        PDK.modToOut(x=subLen, modlen=modLen, isolation=False).put()
        y += wgSpace*2
    
    # Bend radius, 1.5 µm WG
    for br in [200, 100, 50, 25]:
        PDK.inputLabel(text='1.5 µm WG: %.0f µm BR' % br, straight=straight, xs='wgBend').put(0, y)
        PDK.icBend.strt(300 - PDK.sinebend_minlen(wgSpace, br), arrow=False).put()
        PDK.icBend.sinebend(distance=PDK.sinebend_minlen(wgSpace, br), offset=-wgSpace, arrow=False).put()
        PDK.icBend.sinebend(distance=PDK.sinebend_minlen(wgSpace, br), offset=-wgSpace, arrow=False).put(flip=True)
        PDK.modToOut(xs='wgBend', x=subLen, modlen=modLen, isolation=False).put()
        y += wgSpace*2
    
    # SDT
    shallow = 400
    for n in [1, 2, 3]:
        PDK.inputLabel(text='SDT: %i, Shallow: %i' % (2*n, shallow), straight=0).put(0, y)
        for _ in range(n):
            PDK.segmentSDT().flip().put()
            PDK.icShallow.strt(length=shallow/n).put()
            PDK.segmentSDT().put()
            PDK.icPassive.strt(length=wgOverlap).put()
        PDK.modToOut(x=subLen, modlen=modLen, isolation=False).put()
        y += wgSpace*2
    
    # Ecorner
    for n in [1, 2, 3]:
        PDK.inputLabel(text='Ecorner: %i' % (4*n), straight=400).put(0, y)
        for _ in range(n):
            PDK.componentEcorner().put()
            PDK.icPassive.strt(length=25).put()
            PDK.componentEcorner().flip().put()
            PDK.icPassive.strt(length=25).put()
            PDK.componentEcorner().flip().put()
            PDK.icPassive.strt(length=25).put()
            PDK.componentEcorner().put()
            PDK.icPassive.strt(length=25).put()
        PDK.modToOut(x=subLen, modlen=modLen, isolation=False).put()
        y += wgSpace*2
    
    # Ubend
    for n in [1, 2, 3]:
        PDK.inputLabel(text='Ubend: %i' % (2*n), straight=250 + 150*n).put(0, y)
        for _ in range(n):
            PDK.componentUbend().put()
            PDK.icPassive.strt(length=100).put()
            PDK.componentUbend().flip().put()
            PDK.icPassive.strt(length=100).put()
        PDK.modToOut(x=subLen, modlen=modLen, isolation=False).put()
        y += wgSpace*2
    
    # MMI tree for polarization test
    PDK.inputLabel(text='Passive Tree', ssc=True, straight=500).put(0, y)
    PDK.componentFilter().put()
    M = PDK.deviceMMI_tree(N=2, k=2, pitch=20, xs='wgPassive').put()
    O = PDK.straightToOut(x=subLen)
    for i in range(4):
        O.put(M.pin['b'+str(i)])
    y += wgSpace*2



### Assemble die
def die(dieName=dieName):
    # Deduplicate
    if dieName in nazca.cfg.cellnames.keys():
        return nazca.cfg.cellnames[dieName]
    
    with nazca.Cell(instantiate=True, name=dieName) as die:
        ## Utility marks
        # Die outline for lightfield and darkfield
        marks.utility.die(dieSize).put(0,0)
        marks.utility.die(dieSize, layer=1005, grow=100).put(0,0)
        
        # Corner marks
        marks.fCommon.putCorners(PDK.markCorner.remove_layer(['ProtectRidge']), diesize=dieSize, inset=[0,0], flip=True, flop=True)
        
        # MLA150 alignment marks in explicit locations
        for x in [dieSize[0][0]+150, dieSize[1][1]-150]:
            for y in [dieSize[1][1]-150]:
                marks.mla150.AlignHighMag(layer=['ProtectRegrowth', 'ProtectRidge']).put(x,y)
                nazca.netlist.Annotation(layer='Annotation', 
                                        text=( 'Align: %.0f, %.0f' % (x,y) )).put(x,y)
        
        # Stepper2 DFAS marks in explicit locations
        for x in [dieSize[0][0]+350, dieSize[1][1]-350]:
            for y in [dieSize[1][1]-150]:
                marks.stepper2.LocalAlign(layer=['ProtectRidge']).put(x,y)
                nazca.netlist.Annotation(layer='Annotation', 
                                        text=( 'Align: %.0f, %.0f' % (x,y) )).put(x,y)
        
        # Custom reduced-layer resolution tests and verniers in top left
        PDK.cellArray([marks.utility.ResolutionBlockSet(xs_pattern='ProtectRidge'), 
                    marks.utility.ResolutionBlockSet(xs_pattern='EtchRib',  xs_background='ProtectRidge'), 
                    marks.utility.ResolutionBlockSet(xs_pattern='EtchIso',  xs_background='ProtectRidge'), 
                    marks.utility.ResolutionBlockSet(xs_pattern='MetalVia', xs_background='ProtectRidge'), 
                    marks.utility.ResolutionBlockSet(xs_pattern='MetalTop', xs_background='ProtectRidge', 
                                                        res=[0.8, 1.0, 1.2, 1.4, 1.6, 1.8]), 
                    marks.utility.Vernier(xs_center='ProtectRidge', xs_surround='EtchRib'),
                    marks.utility.Vernier(xs_center='ProtectRidge', xs_surround='EtchIso'),
                    marks.utility.Vernier(xs_center='ProtectRidge', xs_surround='MetalVia'),
                    marks.utility.Vernier(xs_center='ProtectRidge', xs_surround='MetalTop'),
                    marks.utility.Vernier(xs_center='ProtectRidge', xs_surround='CleaveLane'),
                    marks.utility.LayerLabels(layers=[1,2,4, 5])],
                    space=20, nx=12).put(dieSize[1][1]/2, dieSize[1][1]-150)
        
        
        ## Place test structures horizontal and vertical with separate labels
        x, y = 0, 0
        for rot in [0, 90]:
            nazca.cp.goto(x, y, rot)
            nazca.cp.shift(0, 50 + 50)
            groupFP.put()
            
            # Labels
            nazca.cp.shift(75, -50)
            marks.utility.label(text=dieName + f' / {str(rot)} deg', height=50, layer=['MetalTop'], 
                                grow=[0], date=True, origin=['lower', 'left']).put()
            x += dieSize[1][0]
        
        
        ## Facet cleave lanes
        # Vertical
        PDK.icCleave.strt(length=dieSize[1][1] + wgOverlap*2, arrow=False).put(wgOverlap/2, -wgOverlap, 90)
        PDK.icCleave.strt(length=dieSize[1][1] + wgOverlap*2, arrow=False).put(dieSize[1][1] - wgOverlap/2, -wgOverlap, 90)
        PDK.icCleave.strt(length=dieSize[1][1] + wgOverlap*2, arrow=False).put(dieSize[1][0] - wgOverlap/2, -wgOverlap, 90)
        
        # Horizontal
        PDK.icCleave.strt(length=dieSize[1][0] + wgOverlap*2, arrow=False).put(-wgOverlap, wgOverlap/2, 0)
        PDK.icCleave.strt(length=dieSize[1][0] + wgOverlap*2, arrow=False).put(-wgOverlap, dieSize[1][1]-wgOverlap/2, 0)
    return die


## Generate GDS
def makeGDS(dieName=dieName, die=die):
    from datetime import datetime
    dieFile = '.\\Process_Tests\\' + dieName + '_' + datetime.now().strftime('%Y%m%d') + '.gds'
    nazca.export_gds(topcells=die, filename=dieFile, clear=False)
    
    # Postprocess with KLayout
    print('Postprocessing with KLayout...')
    from subprocess import call
    from os.path import abspath
    postDRC = [abspath('..\\..\\..\\..\\Code\\Process and Layout\\nazca\\pdk_Fab6\\Fab6_postprocess.lydrc'),
               abspath('..\\..\\..\\..\\Code\\Process and Layout\\nazca\\pdk_Fab6\\Fab6_postprocess_merge.lydrc')]
    for drc in postDRC:
        call(['C:\\Utilities\\EDA\\KLayout\\klayout_app.exe', '-b', '-r', drc, '-rd', 'input='+abspath(dieFile)])

if __name__ == '__main__':
    makeGDS(dieName + ", #000.00", die(dieName + ", #000.00"))

