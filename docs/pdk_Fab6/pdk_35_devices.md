Module pdk_Fab6.pdk_35_devices
==============================
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

Functions
---------

    
`deviceCCDBR(type='VC', order=4, passivelength=60, gratinglength=300, activelength=500, ploc=0.5, l0=1.03, xsCC='wgShallow', arrow=True, instantiate=True)`
:   Corner-cube + DBR laser of specified type.
    
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

    
`deviceDBR(type='VC', order=4, gratinglength=600, rfrac=0.3, activelength=500, ploc=0.5, l0=1.03, arrow=True, instantiate=True)`
:   DBR with specified grating type.
    
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

    
`deviceDFB(type='GC', order=4, length=1500, rfrac=0.3, plen=[1], l0=1.03, d0=3.25, xs2='wgPassive', w2=None, arrow=True, instantiate=True)`
:   DFB with specified grating type.
    
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

    
`deviceMMI(N=2, pitch=45, radius=None, straight=15, dW=None, dL=None, nr=None, xs='wgBend', instantiate=True, arrow=True)`
:   1xN MMI splitter with fan-out to <pitch> pitch.
    
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

    
`deviceMMI22(pitch=45, radius=None, straight=15, dW=None, dL=None, nr=None, xs='wgBend', instantiate=True, arrow=True)`
:   2x2 MMI splitter with fan-in and fan-out to <pitch> pitch.
    
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

    
`deviceMMI_tree(N=2, k=3, pitch=45, radius=None, straight=15, dW=None, dL=None, nr=None, xs='wgBend', instantiate=True, arrow=True)`
:   1xN^k MMI tree; <k> deep with <N> per stage.
    
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

    
`deviceMZM(length=2000, pitch=None, ploc=0.5, dual=True, xsIO='wgPassive', devicefunc=<function deviceModulator>, arrow=True, instantiate=True)`
:   Mach-Zehnder Modulator: 1x2 MMI, modulators, 2x1 MMI.
    
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

    
`deviceModulator(length=1000, isolation=True, ploc=0, xs1='wgPassive', w1=None, xs2=None, w2=None, arrow=True, instantiate=True)`
:   Straight phase modulator, including isolation regions and top metal. Requires electrical connection at either end.
    
    Args:
        length (float): Length of modulator, excluding required isolation
        isolation (bool): Include isolation sections?
        ploc (float): Fractional length along modulator to place electrical pin 'p0'
        xs1, xs2 (str): Cross-section for input and output
        w1, w2 (float): Widths for xs1 and xs2 if nondefault
        arrow (bool): Draw arrows?
    
    Returns:
        Cell: Modulator cell with appropriate pin connections

    
`devicePMArray(length=2000, pitch=45, radius=None, N=2, k=3, dW=None, dL=None, ploc=None, isolation=True, straightin=None, straightout=None, straightsoa=None, soa=None, soacleave=False, soacleavelen=None, outputpitch=None, xs='wgBend', instantiate=True, arrow=True)`
:   1xN phase modulator array.
    
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

    
`deviceSOA(length=500, ploc=0.5, xs1='wgPassive', w1=None, xs2=None, w2=None, arrow=True, instantiate=True)`
:   Straight active region, including top metal and optionally transitions. Requires electrical connection.
    
    Args:
        length (float): Length of active region
        ploc (float): Fractional length along active region to place electrical pin 'p0'
        xs1, xs2 (str): Cross-section for passive lengths; automatically places SDT if not 'wgActive' or 'wgShallow'
        w1, w2 (float): widths for xs1 and xs2 if nondefault
        arrow (bool): Draw arrows?
    
    Returns:
        Cell: Active waveguide cell with appropriate pin connections