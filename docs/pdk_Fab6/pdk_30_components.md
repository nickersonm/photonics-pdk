Module pdk_Fab6.pdk_30_components
=================================
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

Functions
---------

    
`componentCC(angle=44, width=None, xs='wgShallow', instantiate=True)`
:   Corner cube connected to wgPassive.
    
    Args:
        angle (float): angle of closing waveguides
        width (float): waveguide width
        xs (str): cross-section
    
    Returns:
        Cell: properly pinned waveguide cell

    
`componentEcorner(Lbend=2.5, dx=0.48, dy=-0.6, width=2.0, instantiate=True)`
:   Compact elliptical 90° corner.
    
    Args:
        Lbend (float): scale size
        dx (float): x adjustment
        dy (float): y adjustment
        width (float): waveguide width
    
    Returns:
        Cell: properly pinned waveguide cell

    
`componentFan(N=2, pitchin=4, pitchout=127, straightin=0, straightout=0, xs='wgBend', width=None, radius=None, arrow=True, instantiate=True)`
:   Fan in/out from pitchin to pitchout waveguide spacing.
    
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

    
`componentFilter(dW=None, dL=None, xs='wgPassive', arrow=True, instantiate=True)`
:   1x1 MMI filter: one input, one output, rejects higher order modes.
    
    Args:
        dW (float): Modify width
        dL (float): Modify length
        xs (str): cross-section
    
    Returns:
        Cell: Filter cell with pins a0, b0

    
`componentMMI(N=2, dwMMI=None, dlMMI=None, wMMI=None, lMMI=None, nr=None, l0=1.03, taperIn=5, taperOut=5, xs='wgBend', arrow=True, instantiate=True)`
:   1xN MMI: one input, N outputs.
    
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

    
`componentMMI22(dwMMI=None, dlMMI=None, wMMI=None, lMMI=None, nr=None, l0=1.03, taperIn=5, taperOut=5, xs='wgBend', arrow=True, instantiate=True)`
:   2x2 MMI: 2 inputs, 2 outputs.
    
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

    
`componentSSC(width=5, length=80, xs='wgPassive', straight=30.0, angle=0, instantiate=True)`
:   SSC to couple 2 µm MFD free-space gaussian mode to ridge WG.
    
    Args:
        width (float): FS side width
        length (float): Taper length
        straight (float): Straight section of full-width taper for cleaving
        angle (float): Angled facet
        xs (str): output cross-section
    
    Returns:
        Cell: properly pinned waveguide cell

    
`componentUbend(Lin=7, Lout=5, dx=0, dy=-0.32, width=2.0, instantiate=True)`
:   Compact U-bend reversing direction of wgPassive in 2*Lout horizontal space.
    
    Args:
        Lin (float): input length
        Lout (float): output / side length
        dx (float): x adjustment
        dy (float): y adjustment
        width (float): waveguide width
    
    Returns:
        Cell: properly pinned waveguide cell

    
`segmentGCDFB(order=4, w=None, nA=3.442, nB=3.4402, l0=1.03, instantiate=True)`
:   Gain-coupled active region grating cell; no vias or contacts.

    
`segmentInlinePad(length=200, height=125, instantiate=True)`
:   Short inline pad

    
`segmentIsolation(xs1='wgPassive', xs2='wgPassive', w1=None, w2=None, length=40, instantiate=True)`
:   Electrical isolation segment.

    
`segmentLCDBR(order=4, wA=34, wB=6, nA=3.4402, nB=3.4383, l0=1.03, neff=<function <lambda>>, instantiate=True)`
:   Laterally-coupled shallow-rib DBR unit cell.

    
`segmentPad(width=150, length=200, taper=20, instantiate=True)`
:   Electrical bond pad.

    
`segmentSDT(wTR0=None, wTR1=6.5, wShallow=3, wDeep=2, lTR1=6, lTR2=90, lTRs=30, xs2='wgPassive', w2=None, instantiate=True)`
:   Taper from shallow to deep.

    
`segmentVCDBR(order=4, w=None, nA=3.4241, nB=3.4157, l0=1.03, instantiate=True)`
:   Vertically-coupled deep-ridge DBR unit cell

    
`trOpenP(length=None, width=3.0, height=None, contactWidth=30, angle=60, instantiate=True)`
:   Transition to open p-metal.

    
`trPassive2ActiveDeep(length=None, xs1='wgPassive', w1=None, w2=None, contactWidth=30, angle=60, instantiate=True)`
:   Transition from passive deep to active deep.

    
`trPassive2Iso(length=None, xs1='wgPassive', w1=None, w2=None, instantiate=True)`
:   Transition from deep to isolation.

    
`trPassive2Mod(length=None, xs1='wgPassive', w1=None, w2=None, contactWidth=30, angle=60, instantiate=True)`
:   Transition from deep to modulator.

    
`trTaper(length=None, xs1='wgBend', xs2='wgPassive', w1=None, w2=None, instantiate=True)`
:   Generic passive taper.