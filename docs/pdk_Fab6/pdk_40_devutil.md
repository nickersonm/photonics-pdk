Module pdk_Fab6.pdk_40_devutil
==============================
Device utilities on "Fab6" process.
(c) Michael Nickerson 2023

    Device Utilities:
        inputLabel(): Label plus SSC input
        straightToOut(): Straight waveguide until reaching given length, with output termination
        modToOut(): Connect the current position to the output x-coordinate, with a modulator before termination
        soaToOut(): Connect the current position to the output x-coordinate, with an SOA before termination

Functions
---------

    
`inputLabel(text=None, xs='wgPassive', ssc=True, wSSC=None, straight=73, textheight=35, textlayer='MetalTop', textgrow=0, textx=73, texty=23, angle=0, anglecomp=False)`
:   Label plus SSC input with optional straight waveguide.
    
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

    
`modToOut(x, xs='wgPassive', ssc=True, wSSC=None, node=None, modlen=2000, isolation=True, angle=0, anglecomp=False)`
:   Connect the current position to the output x-coordinate, with a modulator before termination.
    
    Args:
        x (float): X-coordinate to terminate at
        xs (str): Cross-section for waveguide
        wSSC (float): Non-default SSC width
        node (nazca.Node): node/pin to use as initial position
        modlen (float): Modulator length
        angle (float): Angled output facet
        anglecomp (bool): Compensate angled output facet with waveguide rotation

    
`soaToOut(x, xs='wgPassive', ssc=True, wSSC=None, node=None, soalen=500, angle=0, anglecomp=False)`
:   Connect the current position to the output x-coordinate, with an SOA before termination.
    
    Args:
        x (float): X-coordinate to terminate at
        xs (str): Cross-section for waveguide
        wSSC (float): Non-default SSC width
        node (nazca.Node): node/pin to use as initial position
        soalen (float): Active region length
        angle (float): Angled output facet
        anglecomp (bool): Compensate angled output facet with waveguide rotation

    
`straightToOut(x, xs='wgPassive', ssc=True, wSSC=None, node=None, angle=0, anglecomp=False)`
:   Straight waveguide until reaching given x-coordinate, with output termination.
    
    Args:
        x (float): X-coordinate to terminate at
        xs (str): Cross-section for waveguide
        wSSC (float): Non-default SSC width
        node (nazca.Node): node/pin to use as initial position
        angle (float): Angled output facet
        anglecomp (bool): Compensate angled output facet with waveguide rotation