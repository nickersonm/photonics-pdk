#!/usr/bin/env python3
# 
# Provides simple components for Fab6 process
# 
# @Authors: Michael Nickerson
# @email: nickersonm@ece.ucsb.edu
# 2023(c)
# 

"""
Simple components for 'Fab6' process.
(c) Michael Nickerson 2023

    Marks:
        markVernier: vernier marks for all layers
        markResolution: resolution tests for all layers
        markDEKTAK: DEKTAK tests for all layers
        markTLM: concentric TLM for P-metal and N-metal
        markCorner: 100x10x100 µm square bracket to mark corners
    Interconnects defined for each cross-section
    Transitions:
        trOpenP(): open a via and p-metal
        trTaper(): taper from xs1 to xs2
        trPassive2Iso(): passive deep WG to isolation WG
        trPassive2Mod(): passive deep WG to modulation WG
        trShallow2Active: shallow WG to modulation WG
    Segments:
        segmentIsolation(): 50 µm of isolation WG with 5 µm passive WG connection on either side
        segmentSDT(): shallow -> deep taper
        segmentGCDFB(): gain-coupled active region grating
        segmentLCDFB(): laterally-coupled active region grating
        segmentVCDBR(): vertically-coupled DBR grating
        segmentLCDBR(): laterally-coupled DBR grating
        segmentPad(): 150x200 µm bond pad
        segmentInlinePad(): 200x125 µm inline bond pad
    Components:
        componentSSC(): spot size converter
        componentMMI(): 1xN MMI
        componentMMI22(): 2x2 MMI
        componentFilter(): 1x1 MMI for mode filtering
        componentUbend(): compact U-bend
        componentEcorner(): compact elliptical 90° corner
        componentCC(): corner cube
        componentFan(): fan-in/fan-out structure
"""

from math import pi, ceil, tan, sin
from numpy import linspace

import nazca

import marks_NanoFab as marks

from .pdk_20_technology import *
from .pdk_05_functions import sinebend_minlen, cellArray, taperLen


### Settings
# Enable DRC
nazca.cfg.group_connect = True
nazca.pin2pin_drc_on()
nazca.pin2pin_drc_raise(num=0)
nazca.cfg.drc_raise = True

# Assorted other settings
version = {
    'cellname': 'auto',
    'version': '1',
    'owner': 'Michael Nickerson'
}

pinstyle = {
    'size': 1.5,
    'stub_length': 0    # Don't want layers drawn for pins
}
nazca.add_pinstyle(name='default', styledict=pinstyle)

pinstyle = {
    'shape': 'inv_pointer',
    'size': 3.0,
    'stub_length': 0.0  # Don't want layers drawn for pins
}
nazca.add_pinstyle(name='metal', styledict=pinstyle)
metalTop.pinstyle = 'metal'
metalP.pinstyle = 'metal'



### Test structures and alignment marks
# Verniers
with nazca.Cell(name='markVernier') as markVernier:
    """Vernier array for all relevant layers"""
    cellArray([marks.utility.Vernier(xs_center='ProtectRidge', xs_surround='ProtectRegrowth'),
               marks.utility.Vernier(xs_center='ProtectRidge', xs_surround='EtchRib'),
               marks.utility.Vernier(xs_center='ProtectRidge', xs_surround='EtchIso'),
               marks.utility.Vernier(xs_center='ProtectRidge', xs_surround='MetalVia'),
               marks.utility.Vernier(xs_center='ProtectRidge', xs_surround='MetalTop')],
              space=20, nx=10).\
                  put(0, 0, 'org')

# Resolution test, all layers
with nazca.Cell(name='markResolution') as markResolution:
    """Resolution test array for all relevant layers"""
    cellArray([marks.utility.ResolutionBlockSet(xs_pattern='ProtectRegrowth', xs_background=1001),
               marks.utility.ResolutionBlockSet(xs_pattern='ProtectRidge',    xs_background=1001),
               marks.utility.ResolutionBlockSet(xs_pattern='EtchRib',         xs_background='ProtectRidge'),
               marks.utility.ResolutionBlockSet(xs_pattern='MetalVia',        xs_background='ProtectRidge',
                                                res=[0.8, 1.0, 1.2, 1.4, 1.6, 1.8]),
               marks.utility.ResolutionBlockSet(xs_pattern='MetalTop',        xs_background='ProtectRidge',
                                                res=[0.8, 1.0, 1.2, 1.5, 2.0, 2.5])],
              space=20, nx=10).\
                  put(0, 0, 'org')

# DEKTAK tests
with nazca.Cell(name='markDEKTAK') as markDEKTAK:
    """DEKTAK marks for all relevant layers"""
    cellArray([marks.utility.DEKTAK_box(xs_pad='ProtectRegrowth', xs_background=1001),
               marks.utility.DEKTAK_box(xs_pad='ProtectRidge', xs_background=1001),
               marks.utility.DEKTAK_box(xs_pad='EtchRib',   xs_background='ProtectRidge'),
               marks.utility.DEKTAK_box(xs_pad='EtchIso',   xs_background='ProtectRidge'),
               marks.utility.DEKTAK_box(xs_pad='MetalVia',  xs_background='ProtectRidge'),
               marks.utility.DEKTAK_box(xs_pad='MetalTop',  xs_background='ProtectRidge')],
              space=20, nx=10).\
        put(0, 0, 'org')

# TLM measurements
with nazca.Cell(name='markTLM') as markTLM:
    """TLM tests"""
    marks.utility.cellShift(marks.utility.ConcentricTLM(xs_pad=['MetalTop', 'MetalVia'], 
                                                        grow_pad=[2, 0], 
                                                        xs_background='ProtectRidge',
                                                        grow_background=10), 
                            ['center', 'right']).put(-20, 0)
    marks.utility.cellShift(marks.utility.CircularTLM(xs_pad=['MetalTop', 'MetalVia'], 
                                                      grow_pad=[1, 0], buffer=40, 
                                                      xs_background='ProtectRidge', 
                                                      grow_background=typicalBuffer, 
                                                      r0=[8, 14], ratio=[1, 1.5, 2.0],
                                                      ratio_outer=1.5), 
                            ['center', 'left']).put(20, 0)

# Corner mark
with nazca.Cell(name='dieCorner') as markCorner:
    """Mark for die corner"""
    marks.utility.layerPolygon(
        poly=[(0, 0), (100, 0), (100, 10),
              (10, 10), (10, 100), (0, 100)], 
        layers=['ProtectRegrowth', 'ProtectRidge']).put(0, 0)



### Waveguide components
## Interconnect objects for waveguide sections
icPassive     = nazca.interconnects.Interconnect(xs='wgPassive')
icBend        = nazca.interconnects.Interconnect(xs='wgBend')
icIso         = nazca.interconnects.Interconnect(xs='wgIso')
icShallow     = nazca.interconnects.Interconnect(xs='wgShallow')
icCleave      = nazca.interconnects.Interconnect(xs='xsCleaveLane')
icModulator   = nazca.interconnects.Interconnect(xs='wgModulator')
icActive      = nazca.interconnects.Interconnect(xs='wgActive')
icActiveDeep  = nazca.interconnects.Interconnect(xs='wgActiveDeep')
icMetalP      = nazca.interconnects.Interconnect(xs='metalP')
icMetalTop    = nazca.interconnects.Interconnect(xs='metalTop')



### Transition elements
## Open P-metal
#   If changing layers here, update pdk_20_technology.metalP as well
def trOpenP(length=None, width=wgActive.width, 
            height=None, 
            contactWidth=contactWidth, angle=60, 
            instantiate=True):
    """Transition to open p-metal."""
    # Process inputs
    height = width - 2*viaInsetRib if height is None else height
    length = ceil(height/tan(angle*pi/180) * 10)/10 if length is None else length
    
    # Name and deduplicate
    name = 'trOpenP.'+str([length, height, contactWidth, angle])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(name=name, instantiate=instantiate) as tr:
        # Pins and connection
        nazca.Pin(name='a0', xs='metalP', width=contactWidth, type='dc').put(0, 0, 180)
        nazca.Pin(name='b0', xs='metalP', width=contactWidth, type='dc').put(length, 0, 0)
        nazca.Pin(name='p0', xs='metalP', width=contactWidth, type='dc').put(-metalBuffer + contactWidth/tan(angle*pi/180)/2, 0, 180)
        
        # Geometry
        marks.utility.layerPolygon(
            poly=nazca.geometries.trapezoid(length=length, 
                                            height=height, 
                                            position=8, angle1=angle, angle2=90),
            layers='MetalVia').put(length)
        marks.utility.layerPolygon(
            poly=nazca.geometries.trapezoid(length=metalBuffer + contactWidth/tan(angle*pi/180)/2 + length, 
                                            height=contactWidth, 
                                            position=8, angle1=angle, angle2=90),
            layers='MetalTop').put(length)
        
    return tr


## Generic taper
def trTaper(length=None, xs1='wgBend', xs2='wgPassive', 
            w1=None, w2=None, instantiate=True):
    """Generic passive taper."""
    # Process inputs
    ic1 = nazca.interconnects.Interconnect(xs=xs1)
    ic2 = nazca.interconnects.Interconnect(xs=xs2)
    
    w1 = ic1._getwidth(width=w1)
    w2 = ic2._getwidth(width=w2)
    length = max([nazca.get_xsection(xs1).taper, 
                  nazca.get_xsection(xs2).taper, 
                  taperLen(w1, w2)]) if length is None else length
    
    # Name and deduplicate
    name = f'trTaper.{xs1}_{xs2}.' + str([length, w1, w2])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(name=name, instantiate=instantiate) as tr:
        # Pins and connection
        a0 = nazca.Pin(name='a0', xs=xs1, width=w1, type='opt').put(0, 0, 180)
        b0 = nazca.Pin(name='b0', xs=xs2, width=w2, type='opt').put(length)
        nazca.connect_optical_path(a0, b0, length, sigtype='opt')
        
        # Geometry
        ic1.sinebend_dw(length=length, 
                        width1=w1,
                        width2=w2, arrow=False).put(0)
    
    return tr


## Deep to isolation
#   If changing layers here, update pdk_20_technology.wgIso as well
def trPassive2Iso(length=None, xs1='wgPassive', 
                  w1=None, w2=None, instantiate=True):
    """Transition from deep to isolation."""
    # Process inputs
    ic1 = nazca.interconnects.Interconnect(xs=xs1)
    w1 = ic1._getwidth(width=w1)
    w2 = icIso._getwidth(width=w2)
    
    # Add taper if differently sized input
    length = wgIso.width + 2*wgIso.outset if length is None else length
    tpIn = None if w2 == w1 else trTaper(xs1=xs1, xs2='wgIso', w1=w1, w2=w2)
    tpLen = tpIn.geolen() if tpIn is not None else 0
    stLen = max([0, length - tpLen])
    length = tpLen + stLen
    
    # Name and deduplicate
    name = f'tr_{xs1}2Iso.' + str([length, w1, w2])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(instantiate=False) as tr:
        # Pins and connection
        a0 = nazca.Pin(name='a0', xs=xs1, width=w1, type='opt').put(0, 0, 180)
        b0 = nazca.Pin(name='b0', xs='wgIso', width=w2, type='opt').put(length, 0, 0)
        nazca.connect_optical_path(a0, b0, length, sigtype='opt')
        
        # Geometry
        if tpIn is not None:
            tpIn.put()
        icPassive.strt(length=stLen, width=w2, arrow=False).put()
        marks.utility.layerPolygon(
            poly=nazca.geometries.trapezoid(length=wgIso.width + 2*wgIso.outset, 
                                            height=wgIso.width + 2*wgIso.outset, 
                                            position=8, angle1=45, angle2=90),
            layers='EtchIso').put()
    
    return tr.flatten(name=name, instantiate=instantiate)


## Deep to modulator
#   If changing layers here, update pdk_20_technology.wgModulator as well
def trPassive2Mod(length=None, xs1='wgPassive', 
                  w1=None, w2=None, 
                  contactWidth=contactWidth, angle=60, instantiate=True):
    """Transition from deep to modulator."""
    # Process inputs
    ic1 = nazca.interconnects.Interconnect(xs=xs1)
    w1 = ic1._getwidth(width=w1)
    w2 = icModulator._getwidth(width=w2)
    
    # Add taper if differently sized input
    length = wgModulator.width + 2*viaInsetRidge if length is None else length
    tpIn = None if w2 == w1 else trTaper(xs1=xs1, xs2='wgModulator', w1=w1, w2=w2)
    tpLen = tpIn.geolen() if tpIn is not None else 0
    stLen = max([0, length - tpLen, metalBuffer])
    length = tpLen + stLen
    
    # Name and deduplicate
    name = f'tr_{xs1}2Mod.' + str([length, w1, w2, contactWidth, angle])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(instantiate=False) as tr:
        # Pins and connection
        a0 = nazca.Pin(name='a0', xs=xs1, width=w1, type='opt').put(0, 0, 180)
        b0 = nazca.Pin(name='b0', xs='wgModulator', width=w2, type='opt').put(length, 0, 0)
        nazca.connect_optical_path(a0, b0, length, sigtype='opt')
        
        # Geometry
        if tpIn is not None:
            tpIn.put()
        icPassive.strt(length=stLen, width=w2, arrow=False).put()
        trOpen = trOpenP(width=wgModulator.width, height=wgModulator.width - 2*viaInsetRidge, 
                         contactWidth=contactWidth, angle=angle)
        trOpen.put('b0', flop=True, flip=True, drc=False).raise_pins(['p0'])
        tr.openlen = trOpen.geolen()
    
    return tr.flatten(name=name, instantiate=instantiate)


## Deep to active deep
#   If changing layers here, update pdk_20_technology.wgActiveDeep as well
def trPassive2ActiveDeep(length=None, xs1='wgPassive', 
                         w1=None, w2=None, 
                         contactWidth=contactWidth, angle=60, instantiate=True):
    """Transition from passive deep to active deep."""
    # Process inputs
    ic1 = nazca.interconnects.Interconnect(xs=xs1)
    w1 = ic1._getwidth(width=w1)
    w2 = icActiveDeep._getwidth(width=w2)
    
    # Add taper if differently sized input
    tpIn = None if w2 == w1 else trTaper(xs1=xs1, xs2='wgActiveDeep', w1=w1, w2=w2)
    tpLen = tpIn.geolen() if tpIn is not None else 0
    length = max([wgActiveDeep.width + 2*viaInsetRidge, metalBuffer + tpLen]) if length is None else length
    stLen = max([0, length - tpLen])
    length = tpLen + stLen
    
    # Name and deduplicate
    name = f'tr_{xs1}2ActiveDeep.' + str([length, w1, w2, contactWidth, angle])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(instantiate=False) as tr:
        # Pins and connection
        a0 = nazca.Pin(name='a0', xs=xs1, width=w1, type='opt').put(0, 0, 180)
        b0 = nazca.Pin(name='b0', xs='wgActiveDeep', width=w2, type='opt').put(length, 0, 0)
        nazca.connect_optical_path(a0, b0, length, sigtype='opt')
        
        # Geometry
        if tpIn is not None:
            tpIn.put()
        icActiveDeep.strt(length=stLen, width=w2, arrow=False).remove_layer(['MetalVia', 'MetalTop']).put()
        trOpen = trOpenP(width=wgActiveDeep.width, height=wgActiveDeep.width - 2*viaInsetRidge, 
                         contactWidth=contactWidth, angle=angle)
        trOpen.put('b0', flop=True, flip=True, drc=False).raise_pins(['p0'])
        tr.openlen = trOpen.geolen()
    
    return tr.flatten(name=name, instantiate=instantiate)


## Passive to active, both shallow
#   If changing layers here, update pdk_20_technology.wgActive as well
with nazca.Cell(instantiate=True, name='trShallow2Active') as trShallow2Active:
    """Transition from passive shallow to active shallow."""
    angle = 60
    height = wgActive.pedestal + 2*typicalBuffer
    length = ceil(height/tan(angle*pi/180) * 10)/10
    
    # Pins and connection
    a0 = nazca.Pin(name='a0', xs='wgShallow', type='opt').put(0, 0, 180)
    b0 = nazca.Pin(name='b0', xs='wgActive', type='opt').put(length)
    nazca.connect_optical_path(a0, b0, length, sigtype='opt')
    
    # Geometry
    icShallow.strt(length=length, width=wgActive.width, arrow=False).put(0)
    trOpen = trOpenP(length=length/2 - typicalBuffer, width=wgShallow.width)
    trOpen.put('b0', flip=True, flop=True, drc=False).raise_pins(['p0'])
    trShallow2Active.openlen = trOpen.geolen()
    marks.utility.layerPolygon(
        poly=nazca.geometries.trapezoid(length=length, 
                                        height=height, 
                                        position=8, angle1=angle, angle2=90),
        layers='ProtectRegrowth').put(length)



### Separate waveguide segments
## Electrical isolation segment
def segmentIsolation(xs1='wgPassive', xs2='wgPassive', 
                     w1=None, w2=None, 
                     length=40, instantiate=True):
    """Electrical isolation segment."""
    # Name and deduplicate
    name = f'SegmentIsolation.{xs1}_{xs2}.'+str([length, w1, w2])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(instantiate=False) as segment:
        tr1 = trPassive2Iso(xs1=xs1, w1=w1)
        tr2 = trPassive2Iso(xs1=xs2, w1=w2)
        tr1.put().raise_pins(['a0'])
        icIso.strt(length, arrow=False).put()
        tr2.flip().put().raise_pins(['b0'])
    
    return segment.flatten(name=name, instantiate=instantiate)


## Shallow-deep taper
def segmentSDT(wTR0=None, wTR1=6.5, wShallow=3, wDeep=2, 
               lTR1=6, lTR2=90, lTRs=30, 
               xs2='wgPassive', w2=None, instantiate=True):
    """Taper from shallow to deep."""
    # Process inputs
    ic = nazca.interconnects.Interconnect(xs=xs2)
    w2 = ic._getwidth(None, w2, 'wgBend')
    wDeep = w2 if wDeep < w2 else wDeep    # Increasing output size is OK without a taper
    wTR0 = wTR0 if wTR0 is not None else wgShallow.pedestal
    
    # Name and deduplicate
    name = f'SegmentSDT.{xs2}.'+str([wTR0, wTR1, wShallow, wDeep, w2, lTR1, lTR2, lTRs])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(instantiate=False) as segment:
        # Deep ridge portions
        icPassive.sinebend_dw(length=lTR1, 
                              width1=wTR0,
                              width2=wTR1, arrow=False).put(0)
        icPassive.sinebend_dw(length=lTR2, 
                              width1=wTR1,
                              width2=wDeep, arrow=False).put()
        
        # Add taper if output needs to shrink
        if w2 < wDeep:
            trTaper(xs1='wgPassive', xs2=xs2, w1=wDeep, w2=w2).put()
        
        # IO pins
        nazca.Pin(name='b0', xs=xs2, width=w2, type='opt').put()
        nazca.Pin(name='a0', xs='wgShallow', width=wShallow, type='opt').put(0, 0, 180)
        
        # Shallow rib portions
        icShallow.sinebend_dw(length=lTRs, 
                              width1=wShallow,
                              width2=wDeep, arrow=False).\
                                  remove_layer('ProtectRidge').put(0, drc=False)
        icShallow.strt(width=wDeep, length=lTR1+lTR2-lTRs, arrow=False).\
            remove_layer('ProtectRidge').put(drc=False)
    
    return segment.flatten(name=name, instantiate=instantiate)


## Gain-coupled active region grating cell; no vias or contacts
def segmentGCDFB(order=4, w=None, bias=0.48, 
                 nA=3.4469, nB=3.4402, l0=1.03, instantiate=True):
    """Gain-coupled active region grating cell; no vias or contacts."""
    # Verify inputs
    w = icActive._getwidth(None, w, 'wgPassive')
    nA = 3.4469 if nA is None else nA
    nB = 3.4402 if nB is None else nB
    
    # Name and deduplicate
    name = 'SegmentGCDFB.'+str([order, w, nA, nB, l0])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    # Calculate parameters
    d = (order + 1/2)*l0 * 2/(nA + nB)
    dA = round(bias*d, 5)
    dB = round((1-bias)*d, 5)
    
    with nazca.Cell(instantiate=False) as segment:
        a0 = nazca.Pin(name='a0', xs='wgActive', type='opt').put(0, 0, 180)
        icShallow.strt(length=dA, width=w, arrow=False).put(0, drc=False)
        icActive.strt(length=dB, width=w, arrow=False).remove_layer(['MetalVia', 'MetalTop']).\
            put(dA, drc=False)
        b0 = nazca.Pin(name='b0', xs='wgActive', type='opt').put(dA+dB, 0, 0)
        nazca.connect_optical_path(a0, b0, 0, sigtype='opt')
    
    segment = segment.flatten(name=name, instantiate=instantiate)
    segment.nA = nA
    segment.nB = nB
    
    return segment


## Vertically-coupled deep-ridge DBR unit cell
def segmentVCDBR(order=5, w=None, bias=0.48, 
                 nA=3.4241, nB=3.4157, l0=1.03, instantiate=True):
    """Vertically-coupled deep-ridge DBR unit cell"""
    # Verify inputs
    w = icPassive._getwidth(None, w, 'wgPassive')
    nA = 3.4241 if nA is None else nA
    nB = 3.4157 if nB is None else nB
    
    # Name and deduplicate
    name = 'SegmentLCDBR.'+str([order, w, nA, nB, l0])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    # Calculate parameters
    d = (order + 1/2)*l0 * 2/(nA + nB)
    dA = round(bias*d, 5)
    dB = round((1-bias)*d, 5)
    
    with nazca.Cell(name=name, instantiate=instantiate) as segment:
        icPassive.strt(width=w, length=dA + dB, arrow=False).put(0, drc=False).raise_pins(['a0', 'b0'])
        marks.utility.layerPolygon(
            poly=nazca.geometries.box(length=dA, width=w+typicalBuffer), 
            layers='EtchRib').put(0)
    
    segment.nA = nA
    segment.nB = nB
    
    return segment


## Laterally-coupled shallow-rib DBR unit cell
def segmentLCDBR(order=3, wA=wgShallow.pedestal, wB=6, bias=0.48, 
                 nA=3.4469, nB=3.4464, l0=1.03, 
                 neff=lambda w: -0.0611*(w**-2.65) + 3.447, instantiate=True):
    """Laterally-coupled shallow-rib DBR unit cell."""
    # If neff lambda provided, use it to recalculate n1 and n2
    if callable(neff):
        nA = round(neff(wA), 6)
        nB = round(neff(wB), 6)
    
    # Verify inputs
    nA = 3.4469 if nA is None else nA
    nB = 3.4464 if nB is None else nB
    
    # Name and deduplicate
    name = 'SegmentLCDBR.'+str([order, wA, nA, wB, nB, l0])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    # Calculate parameters
    d = (order + 1/2)*l0 * 2/(nA + nB)
    dA = round(bias*d, 5)
    dB = round((1-bias)*d, 5)
    
    with nazca.Cell(instantiate=False) as segment:
        icShallow.strt(length=dA+dB, arrow=False).\
            remove_layer('ProtectRidge').\
                put(0, drc=False).raise_pins(['a0', 'b0'])
        icPassive.strt(length=dB, width=wB, arrow=False).put(0,  drc=False)
        icPassive.strt(length=dA, width=wA, arrow=False).put(dB, drc=False)
    
    segment = segment.flatten(name=name, instantiate=instantiate)
    segment.nA = nA
    segment.nB = nB
    
    return segment


## Electrical bond pad
def segmentPad(width=150, length=200, taper=traceWidth, instantiate=True):
    """Electrical bond pad."""
    # Name and deduplicate
    name = 'segmentPad.'+str([width, length, taper])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(instantiate=False) as segment:
        icMetalTop.sinebend_dw(length=length/10, 
                               width1=taper,
                               width2=width, arrow=False).put().raise_pins(['a0'])
        icMetalTop.strt(length=length, width=width, arrow=False).put().raise_pins(['b0'])
        
        # Surround by unetched region
        # icPassive.strt(length=l+taper+20, width=w+20, arrow=False).put(-10, drc=False)
        marks.utility.layerPolygon(poly=segment.get_polygons('MetalTop'), 
                                   layers='ProtectRidge', grow=20, jointype='round').put(0)
    
    return segment.flatten(name=name, instantiate=instantiate)


## Short inline pad
def segmentInlinePad(length=200, height=125, instantiate=True):
    """Short inline pad"""
    # Name and deduplicate
    name = 'segmentInlinePad.'+str([length, height])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(instantiate=False) as segment:
        icMetalTop.strt(length=traceWidth, width=length).put().raise_pins(['a0'])
        icMetalTop.taper(length=height-traceWidth, width1=length, width2=round(length*0.75)).put()
    
    return segment.flatten(name=name, instantiate=instantiate)



### Simple Components
## SSC I/O
def componentSSC(width=5, length=80, xs='wgPassive', 
                 straight=cleaveWidth*1.5, angle=0, 
                 instantiate=True):
    """SSC to couple 2 µm MFD free-space gaussian mode to ridge WG.
    
    Args:
        width (float): FS side width
        length (float): Taper length
        straight (float): Straight section of full-width taper for cleaving
        angle (float): Angled facet
        xs (str): output cross-section
    
    Returns:
        Cell: properly pinned waveguide cell
    """
    if width is None:
        width = 5.0
    ic = nazca.interconnects.Interconnect(xs=xs)
    
    # Name and deduplicate
    name = f'SSC.{xs}.'+str([width, length, straight, angle])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(instantiate=False) as component:
        if angle != 0:
            # Put trapezoid for angled facet
            marks.utility.layerPolygon(
                poly=nazca.geometries.trapezoid(length=round(1.1*width*sin(angle*pi/180), 3), 
                                                height=width, 
                                                position=8, angle1=90-angle, angle2=90),
                layers='ProtectRidge').put()
        ic.strt(length=straight, width=width, arrow=False).put().raise_pins(['a0'])
        ic.sinebend_dw(length=length, 
                       width1=width, 
                       width2=ic.width, 
                       arrow=False).put().raise_pins(['b0'])
    
    return component.flatten(name=name, instantiate=instantiate)


## 1xN MMI
def componentMMI(N=2, dwMMI=None, dlMMI=None, wMMI=None, lMMI=None, 
                 nr=None, l0=1.03, taperIn=5, taperOut=5, 
                 xs='wgBend', arrow=True, instantiate=True):
    """1xN MMI: one input, N outputs.
    
    Args:
        N (int): Number of outputs
        dwMMI (float): Modify width
        dlMMI (float): Modify length; 'None' for optimized default for some N values
        wMMI (float): MMI width; 'None' for calculated
        lMMI (float): MMI length; 'None' for calculated
        nr (float): Slab mode effective index; optionally lambda (l0)
        l0 (float): Design wavelength [µm]
        taperIn (float): Input taper length
        taperOut (float): Output taper length
        xs (str): cross-section for input and output
    
    Returns:
        Cell: MMI cell with pins a0, b0..b(N-1)
    """
    # Process inputs
    N = round(N)    # Force to integer
    ic = nazca.interconnects.Interconnect(xs=xs)
    wWG = ic._getwidth(None, ic.width, 'wgBend')
    nr = 3.441 if nr is None else nr
    l0 = 1.03 if l0 is None else l0
    
    # Hardcoded defaults
    spaceWG = wWG
    dTaper = 0.5 if N > 1 else 0
    
    # Calculate MMI parameters
    if N==3 and wMMI is None and lMMI is None and dlMMI is None:
        # Optimized N=3 MMI parameters
        wMMI = wgPassive.width * 4  # Width for N=2
        lMMI = (nr * wMMI**2) / ( 3 * l0 )   # lPi/4
    if wMMI is None:
        wMMI = N * (spaceWG + wWG)
    dwMMI = 0 if dwMMI is None else dwMMI
    wMMI = wMMI + dwMMI
    if callable(nr):    # If nr = nr(l0)
        nr = nr(l0)
    if nr is None:
        nr = 3.441  # Default AR8.2 epitaxy slab mode index at 1.03 µm
    if taperIn is None:
        taperIn = 5
    if taperOut is None:
        taperOut = 5
    if lMMI is None:
        lMMI = (wMMI**2 * nr) / ( N * l0 )
    
    NdlDefaults = {2.0: {1:-8.5, 2:-9.8, 3:-5.7, 4:-11.5}, 1.5: {1:-7.9, 2:-9.0, 3:-4.7, 4:-10}}
    if dlMMI is None and wWG in NdlDefaults.keys() and N in NdlDefaults[wWG].keys():
        dlMMI = NdlDefaults[wWG][N]
    else:
        dlMMI = 0
    lMMI = lMMI + dlMMI
    
    y = [wMMI * (2*i - (N-1)) / (2*N) for i in range(N)]
    
    # Name and deduplicate
    name = f'MMI.{xs}' + '.1x%i.[%.3g, %.3g, %.3g, %.3g, %.3g]' % (N, wMMI, lMMI, taperIn, taperOut, nr/l0)
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(instantiate=False) as component:
        # Input taper
        ic.sinebend_dw(length=taperIn,
                       width1=wWG, 
                       width2=wMMI, 
                       arrow=False).\
            put().raise_pins(['a0'])
        
        # Main MMI section with a small cut corner at the end
        pinend = ic.strt(length=lMMI-dTaper/2, width=wMMI, arrow=False).put().pinout
        if dTaper > 0:
            pinend = ic.sinebend_dw(length=dTaper/2, 
                                    width1=wMMI, 
                                    width2=wMMI-dTaper/2, 
                                    arrow=False).put().pinout
        
        # Generate and place outputs
        with nazca.Cell(instantiate=False) as outTapers:
            nazca.Pin(name='a0', xs=xs, type='opt', width=wMMI-dTaper/2).put(0,0,180)
            for i, yi in enumerate(y):
                pi = ic.sinebend_dw(length=taperOut, 
                                    width1=(wMMI/N - dTaper/2),
                                    width2=wWG, 
                                    arrow=False).put(0, yi)
                nazca.connect_optical_path(pinend, pi.pinin, 0, sigtype='opt')
                pi.raise_pins(namesin=['b0'], namesout=['b'+str(i)])
        outTapers.put(pinend, drc=None).raise_pins(namesin=['b'+str(i) for i in reversed(range(0,N))])
        
        # Draw pins
        if arrow:
            for _, p in component.ic_pins():
                nazca.make_pincell().put(p)
    
    return component.flatten(name=name, instantiate=instantiate)


## 2x2 MMI
def componentMMI22(dwMMI=None, dlMMI=None, wMMI=None, lMMI=None, 
                    nr=None, l0=1.03, taperIn=5, taperOut=5, 
                    xs='wgBend', arrow=True, instantiate=True):
    """2x2 MMI: 2 inputs, 2 outputs.
    
    Args:
        dwMMI (float): Modify width
        dlMMI (float): Modify length; 'None' for optimized default for some N values
        wMMI (float): MMI width; 'None' for calculated
        lMMI (float): MMI length; 'None' for calculated
        nr (float): Slab mode effective index; optionally lambda (l0)
        l0 (float): Design wavelength [µm]
        taperIn (float): Input taper length
        taperOut (float): Output taper length
        xs (str): cross-section for input and output
    
    Returns:
        Cell: MMI cell with pins a0, b0..b(N-1)
    """
    # Process inputs
    ic = nazca.interconnects.Interconnect(xs=xs)
    wWG = ic._getwidth(None, ic.width, 'wgBend')
    nr = 3.441 if nr is None else nr
    l0 = 1.03 if l0 is None else l0
    
    # Hardcoded defaults
    N = 2
    spaceWG = wWG
    dTaper = 0.5
    
    # Calculate MMI parameters
    if wMMI is None:
        wMMI = N * (spaceWG + wWG)
    dwMMI = 0 if dwMMI is None else dwMMI
    wMMI = wMMI + dwMMI
    if callable(nr):    # If nr = nr(l0)
        nr = nr(l0)
    if nr is None:
        nr = 3.441  # Default AR8.2 epitaxy slab mode index at 1.03 µm
    if taperIn is None:
        taperIn = 5
    if taperOut is None:
        taperOut = 5
    if lMMI is None:
        lMMI = 1.5 * (4 * nr * wMMI**2) / ( 3 * l0 )   # 1.5*lPi
    
    NdlDefaults = {2.0: -12.4, 1.5: -10.8}
    if dlMMI is None and wWG in NdlDefaults.keys():
        dlMMI = NdlDefaults[wWG]
    else:
        dlMMI = 0
    lMMI = lMMI + dlMMI
    
    y = [wMMI * (2*i - (N-1)) / (2*N) for i in range(N)]
    
    # Name and deduplicate
    name = f'MMI.{xs}' + '.2x2.[%.3g, %.3g, %.3g, %.3g, %.3g]' % (wMMI, lMMI, taperIn, taperOut, nr/l0)
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(instantiate=False) as component:
        # Generate and place inputs
        for i, yi in enumerate(y):
            ic.sinebend_dw(length=taperIn, 
                           width1=wWG, 
                           width2=(wMMI/N - dTaper/2),
                           arrow=False).put(0, yi).\
                               raise_pins(namesin=['a0'], namesout=['a'+str(i)])
        
        # Main MMI section with a small cut corner at each end
        pinend = ic.sinebend_dw(length=dTaper/2, 
                                width1=wMMI-dTaper/2, 
                                width2=wMMI, 
                                arrow=False).put(taperIn, 0, drc=False).pinout
        pinend = ic.strt(length=lMMI-dTaper, 
                         width=wMMI, 
                         arrow=False).put().pinout
        pinend = ic.sinebend_dw(length=dTaper/2, 
                                width1=wMMI, 
                                width2=wMMI-dTaper/2, 
                                arrow=False).put().pinout
        
        # Generate and place outputs
        with nazca.Cell(instantiate=False) as outTapers:
            for i, yi in enumerate(y):
                pi = ic.sinebend_dw(length=taperOut, 
                                    width1=(wMMI/N - dTaper/2),
                                    width2=wWG, 
                                    arrow=False).put(0, yi)
                nazca.connect_optical_path(pinend, pi.pinin, 0, sigtype='opt')
                pi.raise_pins(namesin=['b0'], namesout=['b'+str(i)])
        outTapers.put(pinend, drc=None).raise_pins(namesin=['b'+str(i) for i in reversed(range(N))])
        
        # Draw pins
        if arrow:
            for _, p in component.ic_pins():
                nazca.make_pincell().put(p)
    
    return component.flatten(name=name, instantiate=instantiate)


## 1x1 MMI filter
def componentFilter(dW=None, dL=None, xs='wgPassive', 
                    arrow=True, instantiate=True):
    """1x1 MMI filter: one input, one output, rejects higher order modes.
    
    Args:
        dW (float): Modify width
        dL (float): Modify length
        xs (str): cross-section
    
    Returns:
        Cell: Filter cell with pins a0, b0
    """
    # Just return an N=1 MMI
    return componentMMI(N=1, dwMMI=dW, dlMMI=dL, xs=xs, 
                        arrow=arrow, instantiate=instantiate)


## Compact U-bend
def componentUbend(Lin=7, Lout=5, dx=0, dy=-0.32, 
                   width=2.0, instantiate=True):
    """Compact U-bend reversing direction of wgPassive in 2*Lout horizontal space.
    
    Args:
        Lin (float): input length
        Lout (float): output / side length
        dx (float): x adjustment
        dy (float): y adjustment
        width (float): waveguide width
    
    Returns:
        Cell: properly pinned waveguide cell
    """
    # Name and deduplicate
    name = 'Ubend.'+str([Lin, Lout, dx, dy, width])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    # Hardcoded defaults
    w = 0.395*width + 0.1303  # Gaussian MFD conversion
    l0 = 1.03   # Calculated for 1.03 µm only
    nr = 3.441  # Effective slab mode index at l0
    
    # Elliptical reflector calculations
    #   https://doi.org/10.1016/j.optcom.2012.08.078
    #   See EllipticalReflector.nb
    ll = (l0/nr)**2
    Rin = Lin*(1+(pi*w**2/Lin)**2/ll)
    Rout = Lout*(1+(pi**2*w**4 + Lin**2*ll - (pi**4*w**8 + 2*(Lin**2 - 2*Lout**2)*pi**2*w**4*ll + Lin**4*ll**2)**0.5)**2 / (4*Lout*pi*2*w**4*ll))
    a = (Rin + Rout)/2
    b = (Rin*Rout/2)**0.5
    RR = (1 + Rout**2/Rin**2)**0.5
    
    # Polygon calculation
    # Centered-coordinate ellipse
    x = linspace(-a, a, 151)
    y = (b*(a**2 - x**2)**0.5 / a).real
    
    # Rotate pi and shift
    x1 = dx + Lin - Rin/2 + x/RR + Rout*y/(Rin*RR)
    y  = dy + Lout - Rout/2 - Rout*x/(Rin*RR) + y/RR
    x = x1
    
    # Select useful portion
    I = [i for i, _ in enumerate(x) if (y[i] >= 0) and (x[i] >= 0) and (y[i] <= Lout + 3*width)]
    x, y = x[I], y[I]
    
    # Assemble polygon
    maxE = [max(x), max(y)]
    polyE = list(map(list, zip(x, y)))
    polyE = [[0, 0], [0, maxE[1]], *polyE, [maxE[0], 0]]
    
    # Define cell
    with nazca.Cell(instantiate=False) as component:
        # Pins
        nazca.Pin(name='a0', xs='wgPassive', type='opt', width=width).put(0,Lout,180)
        nazca.Pin(name='b0', xs='wgPassive', type='opt', width=width).put(0,-Lout,180)
        
        # Polygon geometry
        polyE = marks.utility.layerPolygon(poly=polyE, 
                                           layers=['ProtectRidge', 'CleaveKeepaway'], 
                                           grow=[0, keepawaySize])
        polyE.put(0)
        polyE.put(0, flip=True)
    
    return component.flatten(name=name, instantiate=instantiate)


## Compact elliptical 90° corner
def componentEcorner(Lbend=2.5, dx=0.48, dy=-0.60, 
                     width=2.0, instantiate=True):
    """Compact elliptical 90° corner.
    
    Args:
        Lbend (float): scale size
        dx (float): x adjustment
        dy (float): y adjustment
        width (float): waveguide width
    
    Returns:
        Cell: properly pinned waveguide cell
    """
    # Name and deduplicate
    name = 'Ecorner.'+str([Lbend, dx, dy, width])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    # Hardcoded defaults
    w = 0.395*width + 0.1303  # Gaussian MFD conversion
    l0 = 1.03   # Calculated for 1.03 µm only
    nr = 3.441  # Effective slab mode index at l0
    
    # Elliptical reflector calculations
    #   https://doi.org/10.1016/j.optcom.2012.08.078
    #   See EllipticalReflector.nb
    rt2 = 2**0.5
    ll = (l0/nr)**2
    R = Lbend*(1+(pi*w**2/Lbend)**2/ll)
    
    # Polygon calculation
    # Centered-coordinate ellipse
    x = linspace(-R, R, 251)
    y = ((R**2 - x**2)**0.5 / rt2).real
    
    # Rotate pi and shift
    x1 = dx + Lbend - R/2 + (x + y)/rt2
    y  = dy - R/2 + (y - x)/rt2
    x = x1
    
    # Select useful portion
    I = [i for i, _ in enumerate(x) if (y[i] >= -Lbend) and (x[i] >= 0) and (y[i] <= 2*width)]
    x, y = x[I], y[I]
    
    # Assemble polygon
    maxE = [max(x), max(y)]
    polyE = list(map(list, zip(x, y)))
    polyE = [[0, -Lbend], [0, maxE[1]], *polyE, [maxE[0], -Lbend]]
    
    # Define cell
    with nazca.Cell(name=name, instantiate=instantiate) as component:
        # Pins
        # TODO: check properly positioned via comparison with Lumerical simulations
        nazca.Pin(name='a0', xs='wgPassive', type='opt', width=width).put(0, 0, 180)
        nazca.Pin(name='b0', xs='wgPassive', type='opt', width=width).put(Lbend, -Lbend, -90)
        
        # Polygon geometry
        marks.utility.layerPolygon(
            poly=polyE, layers=['ProtectRidge', 'CleaveKeepaway'], 
            grow=[0, keepawaySize]).put(0, 0)
    
    return component


## Corner cube
def componentCC(angle=44, width=None, xs='wgShallow', instantiate=True):
    """Corner cube connected to wgPassive.
    
    Args:
        angle (float): angle of closing waveguides
        width (float): waveguide width
        xs (str): cross-section
    
    Returns:
        Cell: properly pinned waveguide cell
    """
    # Process inputs
    ic = nazca.interconnects.Interconnect(xs=xs)
    w = ic._getwidth(None, ic.width, 'wgShallow')
    width = marks.fCommon.cellHeight(ic.strt(length=1).\
        remove_layer(['CleaveKeepaway', 'MetalKeepaway'])) if width is None else width
    
    # Name and deduplicate
    name = f'CC.{xs}.'+str([angle, width])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    # Corner cube geometry
    dX = round(0.5*width/tan(angle*pi/180), 4)
    polyCC = [[0 + dX, 0], [0, width/2], [0, -width/2]]
    
    # Define cell as single straight segment
    with nazca.Cell(instantiate=False) as component:
        ic.strt(length=dX, arrow=False).put(0)
        nazca.Pin(name='a0', xs=xs, width=w).put(0, 0, 180)
        nazca.Pin(name='b0', xs=xs, width=w).put(dX, 0, 0)
    
    # Mask 'ProtectRidge' with polyCC
    component = component.polydiff_layer(clip_poly=polyCC, 
                                         mod_layer='ProtectRidge', 
                                         operation='and')
    
    return component.flatten(name=name, instantiate=instantiate)


## Fan-out/in
def componentFan(N=2, pitchin=4, pitchout=127, 
                 straightin=0, straightout=0, 
                 xs='wgBend', width=None, radius=None, 
                 arrow=True, instantiate=True):
    """Fan in/out from pitchin to pitchout waveguide spacing.
    
    Args:
        N (int): Number of channels
        pitchin (float): input pitch
        pitchout (float): output pitch
        straight(in|out) (float): section of straight waveguide at (in|out)put
        xs (str): cross-section to use
        width (float): waveguide width; default xs.width
        radius (float): minimum radius of curvature; default xs.radius
    
    Returns:
        Cell: properly pinned waveguide cell
    """
    # Process inputs
    N = round(N)
    ic = nazca.interconnects.Interconnect(xs=xs)
    radius = ic._getradius(None, radius, 'wgBend')
    width = ic._getwidth(None, width, 'wgBend')
    
    # Name and deduplicate
    name = 'Fan'+str(N) + '.' + str(xs)+'.'+str([pitchin, pitchout, radius, straightin, straightout, width, int(arrow)])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    # Calculate parameters
    y = [(2*i - (N-1)) / (2*N) for i in range(N)]
    length = ceil( sinebend_minlen(dy=(max(y) * N * abs(pitchout-pitchin)), radius=radius) )
    
    with nazca.Cell(instantiate=False) as component:
        for i, y in enumerate(y):
            ic.strt(length=straightin, arrow=False).put(0, y*pitchin*N).\
                raise_pins(namesin=['a0'],
                           namesout=['a'+str(i)])
            ic.sinebend_dw(length=length, 
                           offset=y*N*(pitchout-pitchin),
                           width=width,
                           N=length*1,
                           arrow=False).put()
            ic.strt(length=straightout, arrow=False).put().\
                raise_pins(namesin=['b0'],
                           namesout=['b'+str(i)])
        
        # Draw pins
        if arrow:
            for _, p in component.ic_pins():
                nazca.make_pincell().put(p)
    
    return component.flatten(name=name, instantiate=instantiate)
