#!/usr/bin/env python3
# 
# Provides complete devices on "Fab6" process
# 
# @Authors: Michael Nickerson
# @email: nickersonm@ece.ucsb.edu
# 2023(c)
# 

"""
Complete devices in 'Fab6' process.
(c) Michael Nickerson 2023

    Devices:
        deviceMMI(): 1xN MMI plus fan-out to specified pitch
        deviceMMI22(): 2x2 MMI plus fan-out to specified pitch
        deviceMMI_tree(): 1xN^k MMI tree, <k> deep with 1xN devices
        deviceModulator(): straight phase modulator
        deviceMZM(): MZM with phase modulators in one or both arms
        devicePMArray(): 1xN phase modulator array
        deviceSOA(): SOA of given length with transitions
        deviceDFB(): DFB laser of specified type
        deviceDBR(): DBR laser of specified type
        deviceCCDBR(): Corner-cube + DBR laser of specified type
"""

from math import floor

import nazca

from .pdk_05_functions import pinGeoDiff
from .pdk_20_technology import *
from .pdk_30_components import *



### Simple Dynamic Devices
## MMI with outputs fanned out
def deviceMMI(N=2, pitch=contactWidth + metalBuffer, radius=None, 
              straight=metalBuffer, dW=None, dL=None, 
              nr=None, xs='wgBend', 
              instantiate=True, arrow=True):
    """1xN MMI splitter with fan-out to <pitch> pitch.
    
    Args:
        N (int): MMI outputs, ≥2
        pitch (float): output pitch
        radius (float): minimum radius of curvature; default xs.radius
        straight (float): straight section before MMI
        dW (float): modify MMI width
        dL (float): modify MMI length
        xs (str): cross-section for input and output
        arrow (bool): draw arrows?
    
    Returns:
        Cell: MMI cell with pins a0, b0..b(N-1)
    """
    # Process inputs
    N = max([2, round(N)])
    ic = nazca.interconnects.Interconnect(xs=xs)
    radius = ic._getradius(None, radius, 'wgBend')
    
    # Name and deduplicate
    name = 'DeviceMMI_1x%i.%s.%s' % (N, xs, str([pitch, float(radius), straight, dW, dL, nr, int(arrow)]) )
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(name=name, instantiate=instantiate) as device:
        ic.strt(length=straight, arrow=False).put().raise_pins(['a0'])
        mmi = componentMMI(N=N, dwMMI=dW, dlMMI=dL, arrow=False, nr=nr, xs=xs).put()
        componentFan(N=N, radius=radius, xs=xs, 
                     pitchin=abs(mmi.pin['b0'].y - mmi.pin['b1'].y), 
                     pitchout=pitch, arrow=False).put().\
                         raise_pins(['b'+str(i) for i in reversed(range(N))])
        
        # Draw pins
        if arrow:
            for _, p in device.ic_pins():
                nazca.make_pincell().put(p)
    
    return device


## 2x2 MMI with outputs fanned out
def deviceMMI22(pitch=contactWidth + metalBuffer, radius=None, 
                straight=metalBuffer, dW=None, dL=None, 
                nr=None, xs='wgBend', 
                instantiate=True, arrow=True):
    """2x2 MMI splitter with fan-in and fan-out to <pitch> pitch.
    
    Args:
        pitch (float): output pitch
        radius (float): minimum radius of curvature; default xs.radius
        straight (float): straight section before MMI
        dW (float): modify MMI width
        dL (float): modify MMI length
        xs (str): cross-section for input and output
        arrow (bool): draw arrows?
    
    Returns:
        Cell: MMI cell with pins a0, b0..b(N-1)
    """
    # Process inputs
    N = 2
    ic = nazca.interconnects.Interconnect(xs=xs)
    radius = ic._getradius(None, radius, 'wgBend')
    
    # Name and deduplicate
    name = 'DeviceMMI_2x2.%s.%s' % (xs, str([pitch, float(radius), straight, dW, dL, nr, int(arrow)]))
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(name=name, instantiate=instantiate) as device:
        mmi = componentMMI22(dwMMI=dW, dlMMI=dL, arrow=False, nr=nr, xs=xs).put()
        componentFan(N=N, radius=radius, xs=xs, 
                     pitchin=pitch, 
                     pitchout=abs(mmi.pin['b0'].y - mmi.pin['b1'].y), 
                     arrow=False).put('b0', mmi.pinin).\
                         raise_pins(['a'+str(i) for i in reversed(range(N))])
        componentFan(N=N, radius=radius, xs=xs, 
                     pitchin=abs(mmi.pin['b0'].y - mmi.pin['b1'].y), 
                     pitchout=pitch, arrow=False).put(mmi.pinout).\
                         raise_pins(['b'+str(i) for i in reversed(range(N))])
        
        # Draw pins
        if arrow:
            for _, p in device.ic_pins():
                nazca.make_pincell().put(p)
    
    return device


## 1xN^k MMI tree
def deviceMMI_tree(N=2, k=3, pitch=contactWidth + metalBuffer, radius=None, 
                   straight=metalBuffer, dW=None, dL=None, nr=None, 
                   xs='wgBend', 
                   instantiate=True, arrow=True):
    """1xN^k MMI tree; <k> deep with <N> per stage.
    
    Args:
        N (int): Number of MMI outputs
        pitch (float): final output pitch
        radius (float): minimum radius of curvature; default xs.radius
        straight (float): straight section before MMI
        dW (float): modify MMI width
        dL (float): modify MMI length
        xs (str): cross-section for input and output
        arrow (bool): Draw arrows?
    
    Returns:
        Cell: MMI cell with pins a0, b0..b(N-1)
    """
    # Process inputs
    N = round(N)
    k = round(k)
    radius = nazca.interconnects.Interconnect(xs=xs)._getradius(None, radius, 'wgBend')
    
    # If no tree needed, return deviceMMI
    if k == 1:
        return deviceMMI(N=N, pitch=pitch, radius=radius, 
                         straight=straight, dW=dW, dL=dL, nr=nr, 
                         xs=xs, instantiate=instantiate, arrow=arrow)
    
    # Name and deduplicate
    name = 'DeviceMMI_1x%i^%i.%s' % (N, k, str([pitch, float(radius), straight, dW, dL, nr, int(arrow)]) )
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(name=name, instantiate=instantiate) as device:
        # Place beginning element
        m = deviceMMI(N=N, pitch=(pitch * N**(k-1)), radius=radius, nr=nr, 
                      dW=dW, dL=dL, straight=straight, xs=xs, arrow=False).put()
        
        # Place sub-trees
        for i in reversed(range(N)):
            deviceMMI_tree(N=N, k=(k-1), pitch=pitch, radius=radius, nr=nr, 
                           straight=straight, dW=dW, dL=dL, xs=xs, arrow=False).\
                put(m.pin['b'+str(i)]).\
                    raise_pins(
                        namesin=['b'+str(j) for j in reversed(range(N**(k-1)))],
                        namesout=['b'+str(N**(k-1)*i + j) for j in reversed(range(N**(k-1)))]
                    )
        
        # Draw pins
        if arrow:
            for _, p in device.ic_pins():
                nazca.make_pincell().put(p)
    
    return device


## Phase modulator
def deviceModulator(length=1000, isolation=True, ploc=0, 
                    xs1='wgPassive', w1=None, 
                    xs2=None, w2=None, 
                    arrow=True, instantiate=True):
    """Straight phase modulator, including isolation regions and top metal. Requires electrical connection at either end.
    
    Args:
        length (float): Length of modulator, excluding required isolation
        isolation (bool): Include isolation sections?
        ploc (float): Fractional length along modulator to place electrical pin 'p0'
        xs1, xs2 (str): Cross-section for input and output
        w1, w2 (float): Widths for xs1 and xs2 if nondefault
        arrow (bool): Draw arrows?
    
    Returns:
        Cell: Modulator cell with appropriate pin connections
    """
    # Preprocess inputs
    xs2 = xs1 if xs2 is None else xs2
    w1 = nazca.interconnects.Interconnect(xs=xs1)._getwidth(None, w1, 'wgPassive')
    w2 = nazca.interconnects.Interconnect(xs=xs2)._getwidth(None, w2, 'wgPassive')
    
    # Name and deduplicate
    name = f'deviceModulator.{xs1}_{xs2}.' + str([length, isolation, ploc, w1, w2, int(arrow)])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    # Just a straight modulator waveguide with isolation segments
    with nazca.Cell(instantiate=False) as device:
        if isolation:
            # Input isolation then passive section
            xsModIO = 'wgPassive'
            wIO = wgPassive.width
            segmentIsolation(xs1=xs1, w1=w1, xs2=xsModIO, w2=wIO).put().raise_pins(['a0'])
        else:
            a0 = nazca.Pin(name='a0', xs=xs1, width=w1, type='opt').put(0, 0, 180)
            xsModIO = xs1
            wIO = w1
        
        # Modulator section with transitions
        p0 = trPassive2Mod(xs1=xsModIO, w1=wIO).put().pin['p0']
        icModulator.strt(length, arrow=False).put()
        p1 = trPassive2Mod(xs1=xsModIO, w1=wIO).flip().put().pin['p0']
        
        if isolation:
            # Passive section then output isolation
            segmentIsolation(xs1=xsModIO, w1=wIO, xs2=xs2, w2=w2).put().raise_pins(['b0'])
        else:
            b0 = nazca.Pin(name='b0', xs=xs2, width=w2, type='opt').put()
            nazca.connect_optical_path(a0, b0, pinGeoDiff(a0, b0) )
        
        # Electrical pins
        nazca.Pin(name='p0', xs='metalTop', type='dc').put(p0.x + ploc*(p1.x - p0.x), 
                                                           contactWidth/2 * int(0 < ploc < 1), 
                                                           180 - 90*int(0 < ploc < 1) - 180*int(ploc == 1))
        nazca.Pin(name='p1', xs='metalTop', type='dc').put(p1.x - ploc*(p1.x - p0.x), 
                                                           -contactWidth/2 * int(0 < ploc < 1), 
                                                           -90*int(0 < ploc < 1) + 180*int(ploc == 1))
        
        # Draw pins
        if arrow:
            for _, p in device.ic_pins():
                nazca.make_pincell().put(p)
    
    return device.flatten(name=name, instantiate=instantiate)


## MZM
def deviceMZM(length=2000, pitch=None, ploc=0.5, 
              dual=True, xsIO='wgPassive', 
              devicefunc=deviceModulator, 
              arrow=True, instantiate=True):
    """Mach-Zehnder Modulator: 1x2 MMI, modulators, 2x1 MMI.
    
    Args:
        length (float): Modulator length
        pitch (float): Spacing between arms
        ploc (float): Fractional distance along modulators for additional electrical connection
        dual (bool): Modulators in both arms
        xs (str): cross-section for passive lengths
        devicefunc (function): Device function for alternative MZMs
        arrow (bool): Draw arrows?
    
    Returns:
        Cell: Modulator cell with appropriate pin connections
    """
    # Process inputs
    ic = nazca.interconnects.Interconnect(xs=xsIO)
    pitch = metalBuffer + contactWidth*(int(dual)+1)/2 if pitch is None else pitch
    
    # Name and deduplicate
    name = f'deviceMZM.{xsIO}.{devicefunc.__name__}.' + str([int(dual)+1, length, pitch, ploc, int(arrow)])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(name=name, instantiate=instantiate) as device:
        # Input MMI and fan-out
        M = deviceMMI(N=2, pitch=pitch, xs=xsIO, arrow=False).put()
        M.raise_pins(['a0'])
        
        # Top arm
        pm = devicefunc(ploc=ploc, length=length, xs1=xsIO, xs2=xsIO, arrow=False)
        pm.put(M.pin['b1']).raise_pins(namesin=['p0'])
        
        # Bottom arm
        if dual:
            pm.put(M.pin['b0'], flip=True).raise_pins(namesin=['p0'], namesout=['p1'])
        else:
            ic.strt(length=pm.geolen(), arrow=False).put(M.pin['b0'])
        
        # Output fan-in and MMI
        deviceMMI(N=2, pitch=pitch, xs=xsIO, arrow=False).put('b0', flip=True).\
            raise_pins(namesin=['a0'], namesout=['b0'])
        
        # Draw pins
        if arrow:
            for _, p in device.ic_pins():
                nazca.make_pincell().put(p)
    
    return device



### More complex static devices
## 1xN phase modulator array
def devicePMArray(length=2000, pitch=contactWidth + metalBuffer, radius=None, 
                  N=2, k=3, dW=None, dL=None, ploc=None, isolation=True, 
                  straightin=None, straightout=None, straightsoa=None, 
                  soa=None, soacleave=False, soacleavelen=None, outputpitch=None, 
                  xs='wgBend', instantiate=True, arrow=True):
    """1xN phase modulator array.
    
    Args:
        length (float): modulated length
        pitch (float): phase modulator pitch
        radius (float): minimum radius of curvature; default xs.radius
        N (int): MMI outputs; total outputs = N^k
        k (int): depth of MMI tree
        dW (float): modify MMI width
        dL (float): modify MMI length
        ploc (float or list): distance along modulators to place electrical pad; if list of length N, applies to individual modulators
        straightin (float): include straight segment of this length before PM
        straightout (float): include straight segment of this length after PM
        straightsoa (float): include straight segment of this length after SOA
        soa (float): include SOA of this length after modulator
        soacleave (bool): include cleave space before SOA; fan-in and fan-out if outputpitch not None
        soacleavelen (float): put SOA cleave at this x-location
        outputpitch (float): pitch for output waveguides
        xs (str): cross-section for passive lengths
    
    Returns:
        Cell: Phase modulator cell with pins a0, b0..b(N-1)
    """
    # Process inputs
    N = round(N)
    k = round(k)
    ic = nazca.interconnects.Interconnect(xs=xs)
    radius = ic._getradius(None, radius, 'wgBend')
    ploc = list(nazca.np.resize([0, 1], N**k)) if ploc is None else ploc
    soa = 0 if soa is None else soa
    soacleavelen = 0 if soacleavelen is None else soacleavelen
    straightout = 0 + ceil(int(soa>0)*(traceWidth*2*(N**k) + metalBuffer)/2) if straightout is None else straightout
    outputpitch = outputpitch if outputpitch is not None else 2*ic.width
    straightout += 3*cleaveWidth if (not outputpitch) and soacleave else 0
    straightin = straightout/2 if straightin is None else straightin
    straightsoa = straightout/2 if straightsoa is None else straightsoa
    
    # Name and deduplicate
    name = f'devicePMArray.{xs}.' + str([N, k, length, pitch, radius, dW, dL, ploc, straightout, soa, int(soacleave), soacleavelen, outputpitch, int(arrow)])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(name=name, instantiate=instantiate) as device:
        # Open with MMI tree
        m = deviceMMI_tree(N=N, k=k, pitch=pitch, radius=radius, 
                           straight=10, arrow=False, xs=xs).put()
        m.raise_pins(['a0'])
        
        for i in reversed(range(N**k)):
            # Get specific ploc
            if not hasattr(ploc, '__iter__'):
                loc = ploc
            elif len(ploc) >= i:
                loc = ploc[i]
            else:
                loc = ploc[0]
            
            # Straight segment before PM
            ic.strt(length=straightin, arrow=False).put(m.pin['b'+str(i)])
            
            # Place modulator
            deviceModulator(length=length, ploc=loc, 
                            xs1=xs, xs2=xs, 
                            isolation=isolation, arrow=False).put().\
                raise_pins(namesin=['p0', 'p1'], namesout=['p'+str(i), 'p'+str(N**k + i)])
            
            # Straight segment between PM and SOA
            ic.strt(length=straightout, arrow=False).put()
            
            # SOA
            if soa > 0 and not soacleave:
                deviceSOA(length=soa, xs1=xs, xs2=xs, ploc=loc, arrow=False).\
                    put().raise_pins(namesin=['p0', 'p1'], namesout=['ap'+str(i), 'ap'+str(N**k + i)])
                
                # Straight segment after SOA
                ic.strt(length=straightsoa, arrow=False).put()
            
            # Output pin if no fan-in
            if outputpitch is None:
                ic.strt(length=0, arrow=False).raise_pins(namesin=['b0'], namesout=['b'+str(i)])
        
        # Fan-in
        if outputpitch is not None:
            F = componentFan(N=N**k, pitchin=pitch, pitchout=outputpitch, 
                             radius=radius, xs=xs, arrow=False)
            FP = F.put()
            
            # SOA if soacleave
            if soa > 0 and soacleave:
                # Generate proper straight length
                soacleavelen = max([cleaveWidth*1.5, soacleavelen - FP.pin['b'+str(i)].x])
                
                # Separation for output cleave
                C = [ic.strt(length=soacleavelen + 1.5*cleaveWidth, arrow=False).put(FP.pin['b'+str(i)]) for i in reversed(range(N**k))]
                nazca.Pin(name='c0', xs=None).put(C[0].pinin.x + soacleavelen, 0, 0)    # Cleave-center pin
                
                # Fan-out again
                FP = componentFan(N=N**k, pitchout=pitch, pitchin=outputpitch, 
                                  radius=radius, xs=xs, arrow=False).put()
                
                # Place SOAs
                for i in reversed(range(N**k)):
                    # Get specific ploc
                    if not hasattr(ploc, '__iter__'):
                        loc = ploc
                    elif len(ploc) >= i:
                        loc = ploc[i]
                    else:
                        loc = ploc[0]
                    deviceSOA(length=soa, xs1=xs, xs2=xs, ploc=loc, arrow=False).\
                        put(FP.pin['b'+str(i)]).raise_pins(namesin=['p0', 'p1'], namesout=['ap'+str(i), 'ap'+str(N**k + i)])
                
                # Fan-in again
                FP = F.put()
            
            # Raise the piuns for the output fan
            FP.raise_pins(namesin=['b'+str(i) for i in range(N**k)])
            
        
        # Draw pins
        if arrow:
            for _, p in device.ic_pins():
                nazca.make_pincell().put(p)
    
    return device


## Straight active region
def deviceSOA(length=500, ploc=0.5, 
              xs1='wgPassive', w1=None, 
              xs2=None, w2=None, 
              arrow=True, instantiate=True):
    """Straight active region, including top metal and optionally transitions. Requires electrical connection.
    
    Args:
        length (float): Length of active region
        ploc (float): Fractional length along active region to place electrical pin 'p0'
        xs1, xs2 (str): Cross-section for passive lengths; automatically places SDT if not 'wgActive' or 'wgShallow'
        w1, w2 (float): widths for xs1 and xs2 if nondefault
        arrow (bool): Draw arrows?
    
    Returns:
        Cell: Active waveguide cell with appropriate pin connections
    """
    # Preprocess inputs
    xs1 = 'wgPassive' if xs1 is None else xs1
    xs2 = xs1 if xs2 is None else xs2
    w1 = nazca.interconnects.Interconnect(xs=xs1)._getwidth(None, w1, 'wgPassive')
    w2 = nazca.interconnects.Interconnect(xs=xs2)._getwidth(None, w2, 'wgPassive')
    
    # Name and deduplicate
    name = f'deviceSOA.{xs1}_{xs2}.' + str([length, ploc, w1, w2, int(arrow)])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    # Just a straight active waveguide with optional transition segments
    with nazca.Cell(instantiate=False) as device:
        if xs1 not in ['wgShallow', 'wgActive']:
            # Input transition
            segmentSDT(xs2=xs1, w2=w1).flip().put(0)
        a0 = nazca.Pin(name='a0', xs=xs1, width=w1, type='opt').put(0, 0, 180)
        
        # Active section with transitions
        p0 = trShallow2Active.put().pin['p0']
        icActive.strt(length - 2*trShallow2Active.openlen, arrow=False).put()
        p1 = trShallow2Active.flip().put(flip=True).pin['p0']
        
        if xs2 not in ['wgShallow', 'wgActive']:
            # Output transition
            segmentSDT(xs2=xs2, w2=w2).put()
        b0 = nazca.Pin(name='b0', xs=xs2, width=w2, type='opt').put()
        nazca.connect_optical_path(a0, b0, pinGeoDiff(a0, b0) )
        
        # Electrical pins
        nazca.Pin(name='p0', xs='metalTop', type='dc').put(p0.x + ploc*(p1.x - p0.x), 
                                                           contactWidth/2 * int(0 < ploc < 1), 
                                                           180 - 90*int(0 < ploc < 1) - 180*int(ploc == 1))
        nazca.Pin(name='p1', xs='metalTop', type='dc').put(p1.x - ploc*(p1.x - p0.x), 
                                                           -contactWidth/2 * int(0 < ploc < 1), 
                                                           -90*int(0 < ploc < 1) + 180*int(ploc == 1))
        
        # Draw pins
        if arrow:
            for _, p in device.ic_pins():
                nazca.make_pincell().put(p)
    
    return device.flatten(name=name, instantiate=instantiate)


## DFB laser of specified type
def deviceDFB(type='GC', order=9, 
              length=1500, rfrac=0.3, 
              plen=[1], l0=1.03, d0=3.25, 
              xs2='wgPassive', w2=None, 
              arrow=True, instantiate=True):
    """DFB with specified grating type.
    
    Args:
        type (string): Type of DFB; only 'GC' available
        order (int): Order of grating; extra whole wavelengths per unit cell
        length (float): Maximum length of DFB [µm]
        rfrac (float): L(right)/L(tot) fraction
        plen (float): Relative size of P-electrodes, will be normalized; 10 µm gaps between
        l0 (float): Design wavelength [µm]
        d0 (float): Between-grating length, in wavelengths
        arrow (bool): Draw arrows?
    
    Returns:
        Cell: Waveguide cell with in and out pins
    """
    # Preprocess inputs
    ic2 = nazca.interconnects.Interconnect(xs=xs2)
    w2 = ic2._getwidth(None, w2, 'wgBend')
    
    # DFB unit cell
    match type:
        case 'GC' | 'GCDFB':
            unitCell=segmentGCDFB(order=order, l0=l0)
        case _:
            nazca.logger.error('Unknown DFB type: ' + type)
    w = unitCell.pinin.width
    d = unitCell.geolen()
    d0 = round(d0*l0/(unitCell.nA), 5) if d0>0 else 0
    
    # DFB parameters
    mR = floor(length*rfrac/d)
    lenR = mR*d
    mL = floor((length-lenR)/d)
    lenL = mL*d
    length = lenR + lenL + d0    # Total length
    
    # Contact sizes
    pGap = 10   # Gap between electrodes
    openP = trOpenP(width=w)
    plen = marks.utility.flatten([plen])
    plen = [round(length*p/sum(plen)/(1 + (len(plen)-1)*pGap/length), 2) for p in plen]
    
    # Name and deduplicate
    name = f'deviceDFB.{xs2}.' + type + '.' + str([order, length, d0, rfrac, plen, l0, w2, int(arrow)])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(name=name, instantiate=instantiate) as device:
        # Place unit cell as left array and right array with quarter-wave waveguide in the middle
        a0 = nazca.Pin(name='a0', xs='wgActive', type='opt').put(0, 0, 180)
        unitCell.flip().put(array=[mL, [d, 0], 1, [0, 0]])
        icActive.strt(length=d0, width=w, arrow=False).remove_layer(['MetalTop', 'MetalVia']).put(mL*d)  # Arrays don't update output pin
        unitCell.put(array=[mR, [d, 0], 1, [0, 0]])
        
        # Output transition
        if xs2 not in ['wgShallow', 'wgActive']:
            sdt = segmentSDT(xs2=xs2, w2=w2).put(length)
            sdt.raise_pins(['b0'])
            b0 = sdt.pinout
        else:
            b0 = nazca.Pin(name='b0', xs='wgActive', type='opt').put(length, 0, 0)
        nazca.connect_optical_path(a0, b0, length, sigtype='opt')
        
        # Electrical connections
        with nazca.Cell(instantiate=False) as contacts:
            for i, lP in enumerate(plen):
                # Vias and top metal
                openP.put()
                icMetalP.strt(width=wgActive.width, length=lP - 2*openP.geolen(), arrow=False).put()
                openP.flip().put(flip=True)
                
                # Contact pins
                nazca.Pin(name='p'+str(i), xs='metalTop', type='dc').put(nazca.cp.x() - lP/2, contactWidth/2, 90)
                nazca.Pin(name='p'+str(len(plen)+i), xs='metalTop', type='dc').put(nazca.cp.x() - lP/2, -contactWidth/2, -90)
                
                # Make gap
                nazca.cp.skip(pGap)
        contacts.flatten().put(0).raise_pins(
            marks.fCommon.flatten([['p'+str(i), 'p'+str(i+len(plen))] for i in range(len(plen))])
        )
        
        if arrow:
            for _, p in device.ic_pins():
                nazca.make_pincell().put(p)
    device.lenR = lenR
    device.lenL = lenL
    
    return device


## DBR laser of specified type
def deviceDBR(type='VC', order=9, 
              gratinglength=600, rfrac=0.30, 
              activelength=500, ploc=0.5, l0=1.03,
              arrow=True, instantiate=True):
    """DBR with specified grating type.
    
    Args:
        type (string): Type of DFB, either 'VC' or 'LC'
        order (int): Order of grating; extra whole wavelengths per unit cell
        gratinglength (float): Maximum length of DBR gratings [µm]
        activelength (float): Length of active region [µm]
        rfrac (float): L(right)/L(tot) fraction
        ploc (float): Fractional length along active region to place electrical pin 'p0'
        l0 (float): Design wavelength [µm]
        arrow (bool): Draw arrows?
    
    Returns:
        Cell: Waveguide cell with in and out pins
    """
    # DBR unit cell
    match type:
        case 'VC' | 'VCDBR':
            unitCell=segmentVCDBR(order=order, l0=l0)
            type = 'VC'
        case 'LC' | 'LCDBR':
            unitCell=segmentLCDBR(order=order, l0=l0)
            type = 'LC'
        case _:
            nazca.logger.error('Unknown DBR type: ' + type)
    d = unitCell.geolen()
    
    # DBR parameters
    mR = floor(gratinglength*rfrac/d)
    lenR = mR*d
    mL = floor((gratinglength-lenR)/d)
    lenL = mL*d
    gratinglength = lenR + lenL
    
    # Name and deduplicate
    name = 'deviceDBR.' + type + '.' + str([order, gratinglength, activelength, rfrac, ploc, l0, int(arrow)])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(name=name, instantiate=instantiate) as device:
        # Place unit cell as left array and right array around SOA
        if mL > 0:
            unitCell.flip().put(array=[mL, [d, 0], 1, [0, 0]])
        soa = deviceSOA(length=activelength, ploc=ploc, 
                        xs1='wgPassive' if mL > 0 and type != 'LC' else 'wgShallow',
                        xs2='wgPassive' if mR > 0 and type != 'LC' else 'wgShallow',
                        arrow=False).put(mL*d)
        soa.raise_pins(['p0', 'p1'])
        if mR > 0:
            unitCell.put(array=[mR, [d, 0], 1, [0, 0]])
        
        # Place pins
        a0 = nazca.Pin(name='a0', xs='wgPassive' if mL > 0 and type != 'LC' else 'wgShallow', type='opt').put(0, 0, 180)
        b0 = nazca.Pin(name='b0', xs='wgPassive' if mR > 0 and type != 'LC' else 'wgShallow', type='opt').put(soa.pinout.x + lenR, 0, 0)
        nazca.connect_optical_path(a0, b0, soa.pinout.x + lenR, sigtype='opt')
        
        if arrow:
            for _, p in device.ic_pins():
                nazca.make_pincell().put(p)
    device.lenR = lenR
    device.lenL = lenL
    
    return device


## Corner-cube + DBR laser of specified type
def deviceCCDBR(type='VC', order=9, passivelength=cleaveWidth*3, 
                gratinglength=300, activelength=500, 
                ploc=0.5, l0=1.03, xsCC='wgShallow', 
                arrow=True, instantiate=True):
    """Corner-cube + DBR laser of specified type.
    
    Args:
        type (string): Type of DFB, either 'VC' or 'LC'
        order (int): Order of grating; extra whole wavelengths per unit cell
        passivelength (float): Length of shallow passive region between CC and active region
        gratinglength (float): Maximum length of DBR output grating [µm]
        activelength (float): Length of active region [µm]
        ploc (float): Fractional length along active region to place electrical pin 'p0'
        l0 (float): Design wavelength [µm]
        xsCC (str): cross-section of corner cube reflector
        arrow (bool): Draw arrows?
    
    Returns:
        Cell: Waveguide cell with in and out pins
    """
    # DBR unit cell
    match type:
        case 'VC' | 'VCDBR':
            unitCell=segmentVCDBR(order=order, l0=l0)
        case 'LC' | 'LCDBR':
            unitCell=segmentLCDBR(order=order, l0=l0)
        case None | 'empty' | 'Empty':
            unitCell=icPassive.strt(length=0)
        case _:
            nazca.logger.error('Unknown DBR type: ' + type)
    d = unitCell.geolen()
    
    # DBR parameters
    m = floor(gratinglength/d)
    gratinglength = m*d
    
    # Name and deduplicate
    name = f'deviceCCDBR.{type}.{xsCC}.' + str([order, gratinglength, activelength, ploc, l0, int(arrow)])
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(name=name, instantiate=instantiate) as device:
        # Place unit cell as left array and right array
        componentCC(xs=xsCC).flip().put(0)
        if xsCC in ['wgShallow', 'wgActive']:
            icShallow.strt(length=passivelength, arrow=False).put()
        else:
            nazca.interconnects.Interconnect(xs=xsCC).strt(length=passivelength, arrow=False).put()
        soa = deviceSOA(length=activelength, ploc=ploc, 
                        xs1=xsCC, xs2='wgPassive' if type != 'LC' else 'wgShallow',
                        arrow=False).put()
        soa.raise_pins(['p0', 'p1'])
        unitCell.put(array=[m, [d, 0], 1, [0, 0]])
        
        # Place pins
        a0 = nazca.Pin(name='a0', xs=xsCC, type='opt').put(0, 0, 180)
        b0 = nazca.Pin(name='b0', xs='wgPassive' if type != 'LC' else 'wgShallow', type='opt').put(soa.pinout.x + gratinglength, 0, 0)
        nazca.connect_optical_path(a0, b0, soa.pinout.x + gratinglength, sigtype='opt')
        
        if arrow:
            for _, p in device.ic_pins():
                nazca.make_pincell().put(p)
    device.lenR = gratinglength
    
    return device


