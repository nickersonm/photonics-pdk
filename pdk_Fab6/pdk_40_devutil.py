#!/usr/bin/env python3
# 
# Provides device utilities on "Fab6" process
# 
# @Authors: Michael Nickerson
# @email: nickersonm@ece.ucsb.edu
# 2023(c)
# 

"""
Device utilities on "Fab6" process.
(c) Michael Nickerson 2023

    Device Utilities:
        inputLabel(): Label plus SSC input
        straightToOut(): Straight waveguide until reaching given length, with output termination
        modToOut(): Connect the current position to the output x-coordinate, with a modulator before termination
        soaToOut(): Connect the current position to the output x-coordinate, with an SOA before termination
"""

import nazca

from .pdk_20_technology import *
from .pdk_30_components import *
from .pdk_35_devices import *

from math import asin, sin


### Device Utilities
## Label plus SSC input
def inputLabel(text=None, xs='wgPassive', ssc=True, wSSC=None, straight=keepawaySize+metalBuffer+typicalBuffer*2, 
               textheight=textHeight, textlayer='MetalTop', textgrow=0, 
               textx=keepawaySize+metalBuffer+typicalBuffer*2, texty=metalBuffer + 2*typicalBuffer,
               angle=0, anglecomp=False):
    """Label plus SSC input with optional straight waveguide.
    
    Args:
        text (str): Label of this waveguide
        xs (str): Cross-section for waveguide; if None then no waveguide is created
        ssc (bool): Include SSC?
        wSSC (float): Non-default SSC width
        straight (float): Straight waveguide length
        textheight (float): Text height
        textlayer (layer): Layer of label
        textgrow (float): Grow text by this amount
        textx (float): Text x-offset
        texty (float): Text y-offset
        angle (float): Angled input facet
        anglecomp (bool): Compensate angled output facet with waveguide rotation
    """
    # Process inputs
    ic = nazca.interconnects.Interconnect(xs=xs) if xs is not None else None
    
    if not ssc and xs is not None:
        wSSC = ic._getwidth(None, ic.width, xs='wgPassive')
    if wSSC is None:
        wSSC = 5.0
    
    with nazca.Cell(instantiate=False) as segment:
        if xs is not None:
            if anglecomp and angle != 0:
                # Calculate waveguide rotation to compensate for etched facet output
                rot = round(asin(3.441*sin(angle*pi/180))*180/pi - 2*angle, 5)
            else:
                # Tilt waveguide so cleave hits at correct angle
                rot = angle
            
            componentSSC(width=wSSC, xs=xs, angle=angle, 
                         straight=1.5*cleaveWidth if rot==angle else 0).put(0, 0, -rot).raise_pins(['a0'])
            if angle != 0:
                ic.euler_bend(angle=rot).put()
            ic.strt(length=straight, arrow=False).put().raise_pins(['b0'])
        else:
            nazca.Pin(name='b0', xs=xs, type='opt').put(straight, 0, 0)
        
        if text is not None:
            marks.utility.label(text=text, height=textheight, date=False, grow=textgrow, 
                                layer=textlayer, origin=['lower', 'left']).\
                                    put(textx, texty)
    return segment


## Straight waveguide until reaching given length, with output termination
def straightToOut(x, xs='wgPassive', ssc=True, wSSC=None, 
                  node=None, angle=0, anglecomp=False):
    """Straight waveguide until reaching given x-coordinate, with output termination.
    
    Args:
        x (float): X-coordinate to terminate at
        xs (str): Cross-section for waveguide
        wSSC (float): Non-default SSC width
        node (nazca.Node): node/pin to use as initial position
        angle (float): Angled output facet
        anglecomp (bool): Compensate angled output facet with waveguide rotation
    """
    # Process inputs
    ic = nazca.interconnects.Interconnect(xs=xs)
    x0 = nazca.cp.x() if node is None else (node.x() if callable(node.x) else node.x)
    
    with nazca.Cell(instantiate=False) as segment:
        inputLabel(ssc=ssc, wSSC=wSSC, xs=xs, angle=angle, anglecomp=anglecomp).\
            put(x - x0, 0, 180).raise_pins(namesin=['a0'], namesout=['b0'])
        ic.strt(length=nazca.cp.x(), arrow=False).\
            put().raise_pins(namesin=['b0'], namesout=['a0'])
    
    return segment


## Place a modulator then take a waveguide to the output
def modToOut(x, xs='wgPassive', ssc=True, wSSC=None, node=None, 
             modlen=2000, isolation=True, angle=0, anglecomp=False):
    """Connect the current position to the output x-coordinate, with a modulator before termination.
    
    Args:
        x (float): X-coordinate to terminate at
        xs (str): Cross-section for waveguide
        wSSC (float): Non-default SSC width
        node (nazca.Node): node/pin to use as initial position
        modlen (float): Modulator length
        angle (float): Angled output facet
        anglecomp (bool): Compensate angled output facet with waveguide rotation
    """
    # Process inputs
    ic = nazca.interconnects.Interconnect(xs=xs)
    x0 = nazca.cp.x() if node is None else (node.x() if callable(node.x) else node.x)
    
    with nazca.Cell(instantiate=False) as segment:
        inputLabel(ssc=ssc, wSSC=wSSC, straight=cleaveWidth, xs=xs, angle=angle, anglecomp=anglecomp).\
            put(x - x0, 0, 180).raise_pins(namesin=['a0'], namesout=['b0'])
        deviceModulator(length=modlen, xs1=xs, isolation=isolation, arrow=False).put(flip=True).raise_pins(['p0'])
        if nazca.cp.get_xy()[0] > 0:
            ic.strt(length=nazca.cp.x(), arrow=False).\
                put().raise_pins(namesin=['b0'], namesout=['a0'])
    
    return segment


## Connect the current position to the output x-coordinate, with an SOA before termination
def soaToOut(x, xs='wgPassive', ssc=True, wSSC=None, node=None, 
             soalen=500, angle=0, anglecomp=False):
    """Connect the current position to the output x-coordinate, with an SOA before termination.
    
    Args:
        x (float): X-coordinate to terminate at
        xs (str): Cross-section for waveguide
        wSSC (float): Non-default SSC width
        node (nazca.Node): node/pin to use as initial position
        soalen (float): Active region length
        angle (float): Angled output facet
        anglecomp (bool): Compensate angled output facet with waveguide rotation
    """
    # Process inputs
    ic = nazca.interconnects.Interconnect(xs=xs)
    x0 = nazca.cp.x() if node is None else (node.x() if callable(node.x) else node.x)
    
    with nazca.Cell(instantiate=False) as segment:
        inputLabel(ssc=ssc, wSSC=wSSC, straight=textHeight, angle=angle, anglecomp=anglecomp).\
            put(x - x0, 0, 180).raise_pins(namesin=['a0'], namesout=['b0'])
        deviceSOA(length=soalen, xs1=xs, xs2=xs, arrow=False).put(flip=True).raise_pins(['p0'])
        if nazca.cp.get_xy()[0] > 0:
            ic.strt(length=nazca.cp.x(), arrow=False).\
                put().raise_pins(namesin=['b0'], namesout=['a0'])
    
    return segment

