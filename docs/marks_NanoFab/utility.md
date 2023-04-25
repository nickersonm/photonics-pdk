Module marks_NanoFab.utility
============================
Utility / diagnostic markings.
(c) Michael Nickerson 2023

Functions
---------

    
`CircularTLM(xs_pad, xs_background=None, grow_pad=None, grow_background=10, r0=[5, 6.5, 9], ratio=[1, 1.5, 2.5], buffer=40, ratio_outer=1.25)`
:   Array of circular TLM patterns as in Pan 2013 (<https://doi.org/10/ghjtgf>).
    
    Args:
        xs_pad (xs or layer): Pads defined in this cross-section [str] or layer(s) [int or list(int)]
        xs_background (xs or layer, optional): Background defined in this cross-section or layer
        grow_* (int or list(int), optional): Grow associated geometry accordingly
        r0 (int, optional): Initial central dot size; repeats along y
        ratio (list, optional): Ratio between successive radii along x. Should have 3 elements.
        buffer (int, optional): Buffer between patterns and between inner and outer contact edges
        ratio_outer (int, optional): Ratio of outer ring electrode to inner electrode
    
    Returns:
        Cell: Left-pad-origin CircularTLM pad array

    
`ConcentricTLM(xs_pad, xs_background=None, grow_pad=None, grow_background=20, r0=15, ratio=[1.65, 2.74, 4.34, 5.45])`
:   Concentric circular-pad TLM pattern as in Reeves 1980 (<https://doi.org/10/dqkw2f>).
    
    Args:
        xs_pad (xs or layer): Pads defined in this cross-section [str] or layer(s) [int or list(int)]
        xs_background (xs or layer, optional): Background defined in this cross-section or layer
        grow_* (int or list(int), optional): Grow associated geometry accordingly
        r0 (int, optional): Center dot size.
        ratio (list, optional): Ratio between successive radii. Defaults to Reeves 1980 values. Must have 4 elements.
    
    Returns:
        Cell: Left-pad-origin ConcentricTLM pad array

    
`DEKTAK(xs_pad, xs_background=None, grow_pad=None, grow_background=10)`
:   Three 20x80 µm pads separated by 20 µm for DEKTAK depth measurement.
    
    Args:
        xs_pad (xs or layer): Pads defined in this cross-section [str] or layer(s) [int or list(int)].
        xs_background (xs or layer, optional): Background defined in this cross-section or layer.
        grow_* (int or list(int), optional): Grow associated geometry accordingly
    
    Returns:
        Cell: Centered-origin DEKTAK measurement mark

    
`DEKTAK_box(xs_pad, xs_background=None, grow_pad=None, grow_background=10)`
:   20 µm thick box/frame across 80x80 µm for DEKTAK depth measurement.
    
    Args:
        xs_pad (xs or layer): Pads defined in this cross-section [str] or layer(s) [int or list(int)].
        xs_background (xs or layer, optional): Background defined in this cross-section or layer.
        grow_* (int or list(int), optional): Grow associated geometry accordingly
    
    Returns:
        Cell: Centered-origin DEKTAK measurement mark

    
`LayerLabels(layers=None, size=25, labels=[])`
:   Generates one or more labels for passed layers, in the form of '[#] layername'.
    
    Args:
        layers (int, str, or list, optional): Layer(s) to generate labels for; if None, will use all layers
        size (int, optional): Text size (height). Defaults to 25.
        labels (list of str, optional): Labels to use instead of nazca-defined layer name
    
    Returns:
        Cell: Left-centered-origin cell of layer labels

    
`ResolutionBlock(xs_pattern, xs_background=None, res=0.5, box_size=15, grow_pattern=None, grow_background=5)`
:   Resolution test pattern; all quadrants visible for test resolution, 
            diagonals for lower resolution, one side for uneven thresholding.
    
    Args:
        xs_pattern (xs or layer): Cross-section [str] or layer(s) [int or list(int)] for test pattern.
        xs_background (xs or layer, optional): Cross-section or layer for background box.
        res (float): Resolution to test in µm.
        box_size (float): Test quadrant size in µm.
        grow_* (int or list(int), optional): Grow associated geometry accordingly.
    
    Returns:
        Cell: Centered-origin resolution test pattern

    
`ResolutionBlockSet(xs_pattern, res=[0.6, 0.7, 0.8, 0.9, 1.0, 1.1], box_size=15, with_label=True, label_size=20, xs_background=None, grow_pattern=None, grow_background=5)`
:   Set of resolution block test patterns at different resolutions, with labels, in two rows.
    
    Args:
        xs_pattern (xs or layer): Cross-section [str] or layer(s) [int or list(int)] for test pattern.
        res (float): Resolution(s) to test in µm.
        box_size (float): Test quadrant size in µm.
        with_label (bool): Include resolution label
        label_size (float): Height of label text
        xs_background (xs or layer, optional): Cross-section or layer for background box.
        grow_* (int or list(int), optional): Grow associated geometry accordingly.
    
    Returns:
        Cell: Centered-origin set of resolution test patterns in two rows

    
`SquareTLM(xs_pad, xs_background=None, grow_pad=None, grow_background=10, padsize=100, spacing=[3, 4, 6, 10, 15, 20, 25])`
:   Square-pad TLM pattern.
    
    Args:
        xs_pad (xs or layer): Pads defined in this cross-section [str] or layer(s) [int or list(int)]
        xs_background (xs or layer, optional): Background defined in this cross-section or layer
        grow_* (int or list(int), optional): Grow associated geometry accordingly
        padsize (int, optional): Length of pad side. Defaults to 100.
        spacing (list, optional): Space between pads. Defaults to [3, 4, 6, 10, 15, 20, 25].
    
    Returns:
        Cell: Left-pad-origin TLM pad array

    
`Vernier(xs_center, xs_surround, xs_background=None, grow_background=10)`
:   Vernier marks for layer alignment verification.
    
    Args:
        xs_center (xs or layer): Cross-section [str] or layer(s) for center mark
        xs_surround (xs or layer): Cross-section [str] or layer(s) for surrounding mark
        xs_background (xs or layer, optional): Bounding background defined in this cross-section or layer
        grow_* (int or list(int), optional): Grow associated geometry accordingly
    
    Returns:
        Cell: Centered-origin vernier mark

    
`die(size, layer=1001, grow=None)`
:   Generates die outline.
    
    Args:
        size (list<float>): size dxy or [dx,dy] for centered die, or lower-left and upper-right corner [[x0, y0], [x1, y1]]
        layer (int, optional): Outline layer; likely 'comment' layer. Defaults to 1001.
    
    Returns:
        Cell: Centered-origin outline

    
`label(text, height=50, layer=1, grow=0, origin=['lower', 'center'], date=False)`
: