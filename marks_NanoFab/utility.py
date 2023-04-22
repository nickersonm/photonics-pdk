#!/usr/bin/env python3
# 
# Provides standard UCSB NanoFab tool markings
#   This file for typical utility (diagnostic) marks
# 
# @Authors: Michael Nickerson
# @email: nickersonm@ece.ucsb.edu
# 2020(c)
# 
# TODO later: merge various outputs with fCommon.polygonMerge

"""
Utility / diagnostic markings
"""

import nazca
from numpy.lib.function_base import flip
from math import floor
from . import pathGDS
from .fCommon import *
from numpy import mod
from datetime import datetime

# Utility marks
# Label, optionally with date
def label(text, height=50, layer=1, grow=0, origin=['lower', 'center'], date=False):
    # Name
    name = 'text_'+text+'.'+str([height, layer, grow, origin, date])
    
    c = nazca.text(text=text + (' (' + datetime.now().strftime('%Y%m%d') + ')' if date else ''), 
                   align='lc', height=height, layer=1000)
    c = cellShift(layerPolygon(layers=layer,
                               poly=cellPolygons(c, 1000),
                               grow=grow), origin)
    return c.flatten(name=name, instantiate=True)

# Simple box to define a die
def die(size, layer=1001, grow=None):
    """
    die Generates die outline

    Args:
        size (list<float>): size dxy or [dx,dy] for centered die, or lower-left and upper-right corner [[x0, y0], [x1, y1]]
        layer (int, optional): Outline layer; likely 'comment' layer. Defaults to 1001.
    
    Returns:
        Cell: Centered-origin outline
    """
    # Process input
    if len(size) == 1:   # Square input
        size = [size, size]
    if isinstance(size[0], float) or isinstance(size[0], int):    # Centered rectangle
        size = [[-size[0]/2, -size[1]/2], [size[0]/2, size[1]/2]]
    
    name = 'die.'+str(layer)+'.'+str(size)
    if name in nazca.cfg.cellnames.keys():
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(name=name) as mark:
        layerPolygon(layers=layer, poly=[size[0], 
                                         [size[0][0], size[1][1]],
                                         size[1],
                                         [size[1][0], size[0][1]]], grow=grow).put(0)
    
    return mark


# Typical style DEKTAK triple-ridge
def DEKTAK(xs_pad, xs_background=None, grow_pad=None, grow_background=10):
    """DEKTAK:
        Three 20x80 µm pads separated by 20 µm for DEKTAK depth measurement

    Args:
        xs_pad (xs or layer): Pads defined in this cross-section [str] or layer(s) [int or list(int)].
        xs_background (xs or layer, optional): Background defined in this cross-section or layer.
        grow_* (int or list(int), optional): Grow associated geometry accordingly

    Returns:
        Cell: Centered-origin DEKTAK measurement mark
    """
    # Name and deduplicate
    name = 'DEKTAK.'+str(layerlist(xs_pad))
    name += '/'+str(layerlist(xs_background)) if xs_background is not None else ''
    if name in nazca.cfg.cellnames.keys():
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(name=name) as mark:
        if (xs_background is not None):
            layerRectangle(layers=xs_background, dx=100, dy=80, grow_layers=grow_background).put(0)
        for x in [-40,0,40]:
            layerRectangle(layers=xs_pad, dx=20, dy=80, grow_layers=grow_pad).put(x)
    return mark


# DEKTAK box
def DEKTAK_box(xs_pad, xs_background=None, grow_pad=None, grow_background=10):
    """DEKTAK_box:
        20 µm thick box/frame across 80x80 µm for DEKTAK depth measurement

    Args:
        xs_pad (xs or layer): Pads defined in this cross-section [str] or layer(s) [int or list(int)].
        xs_background (xs or layer, optional): Background defined in this cross-section or layer.
        grow_* (int or list(int), optional): Grow associated geometry accordingly

    Returns:
        Cell: Centered-origin DEKTAK measurement mark
    """
    # Name and deduplicate
    name = 'DEKTAK_box.'+str(layerlist(xs_pad))
    name += '/'+str(layerlist(xs_background)) if xs_background is not None else ''
    if name in nazca.cfg.cellnames.keys():
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(name=name) as mark:
        if (xs_background is not None):
            layerRectangle(layers=xs_background, dx=100, grow_layers=grow_background).put(0)
        cellShift(layerPolygon(layers=xs_pad, 
                               poly=nazca.geom.frame(sizew=20, sizel=60, sizeh=60), 
                               grow=grow_pad), 'center').put(0,0)
    return mark


# Vernier
def Vernier(xs_center, xs_surround, xs_background=None, grow_background=10):
    """Vernier:
        Vernier marks for layer alignment verification

    Args:
        xs_center (xs or layer): Cross-section [str] or layer(s) for center mark
        xs_surround (xs or layer): Cross-section [str] or layer(s) for surrounding mark
        xs_background (xs or layer, optional): Bounding background defined in this cross-section or layer
        grow_* (int or list(int), optional): Grow associated geometry accordingly

    Returns:
        Cell: Centered-origin vernier mark
    """
    xs_center = flatten(layerlist(xs_center))
    xs_surround = flatten(layerlist(xs_surround))
    
    # Name and deduplicate
    name = 'Vernier'+str(xs_center)+'.'+str(xs_surround)
    name += '/'+str(layerlist(xs_background)) if xs_background is not None else ''
    if name in nazca.cfg.cellnames.keys():
        return nazca.cfg.cellnames[name]
    
    # TODO: make own vernier scale: 0.1 µm/mark short misalignment on inside scale
    # Load GDS
    vernier = nazca.cfg.cellnames['StdVernier'] if 'StdVernier' in nazca.cfg.cellnames.keys() else \
              nazca.load_gds(filename=pathGDS+'Vernier.gds', instantiate=False, 
                             cellname='vernier', newcellname='StdVernier')

    with nazca.Cell(instantiate=False) as mark:
        for l in xs_center:
            vernier.rebuild(flat=True, layermapmode='none', 
                            layermap={1:l}, instantiate=False).put(0)
        for l in xs_surround:
            vernier.rebuild(flat=True, layermapmode='none', 
                            layermap={2:l}, instantiate=False).put(0)
        
        if (xs_background is not None):
            layerPolygon(layers=xs_background, poly=vernier.bbox_polygon.tolist(), grow=grow_background).put(0)
        
    return mark.flatten(name=name, instantiate=True)


# Square TLM pads
def SquareTLM(xs_pad, xs_background=None, grow_pad=None, grow_background=10, padsize=100, spacing=[3, 4, 6, 10, 15, 20, 25]):
    """
    SquareTLM:
        Square-pad TLM pattern

    Args:
        xs_pad (xs or layer): Pads defined in this cross-section [str] or layer(s) [int or list(int)]
        xs_background (xs or layer, optional): Background defined in this cross-section or layer
        grow_* (int or list(int), optional): Grow associated geometry accordingly
        padsize (int, optional): Length of pad side. Defaults to 100.
        spacing (list, optional): Space between pads. Defaults to [3, 4, 6, 10, 15, 20, 25].

    Returns:
        Cell: Left-pad-origin TLM pad array
    """
    # Name and deduplicate
    name = 'TLM.Square.'+str(layerlist(xs_pad))
    name += '/'+str(layerlist(xs_background)) if xs_background is not None else ''
    if name in nazca.cfg.cellnames.keys():
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(name=name) as mark:
        with nazca.Cell(instantiate=False) as pads:
            x0 = 0
            for s in spacing:
                layerRectangle(layers=xs_pad, dx=padsize, grow_layers=grow_pad).put(x0)
                x0 += padsize + s
        if (xs_background is not None):
            layerPolygon(layers=xs_background, poly=pads.bbox_polygon.tolist(), grow=grow_background).put(0)
        pads.put(0)
    
    return mark


# Concentric TLMs [Reeves 1980 https://doi.org/10/dqkw2f ]
def ConcentricTLM(xs_pad, xs_background=None, grow_pad=None, grow_background=20, r0=15, ratio=[1.65, 2.74, 4.34, 5.45]):
    """
    ConcentricTLM:
        Concentric circular-pad TLM pattern as in Reeves 1980 ( https://doi.org/10/dqkw2f )

    Args:
        xs_pad (xs or layer): Pads defined in this cross-section [str] or layer(s) [int or list(int)]
        xs_background (xs or layer, optional): Background defined in this cross-section or layer
        grow_* (int or list(int), optional): Grow associated geometry accordingly
        r0 (int, optional): Center dot size.
        ratio (list, optional): Ratio between successive radii. Defaults to Reeves 1980 values. Must have 4 elements.

    Returns:
        Cell: Left-pad-origin ConcentricTLM pad array
    """
    if len(ratio) < 4:
        nazca.main_logger(msg='"ratio" does not have required 4 elements!', level='error')
    
    # Name and deduplicate
    name = 'TLM.Concentric.'+str(r0)+'.'+str(layerlist(xs_pad))
    name += '/'+str(layerlist(xs_background)) if xs_background is not None else ''
    if name in nazca.cfg.cellnames.keys():
        return nazca.cfg.cellnames[name]
    
    [r1p, r1, r2p, r2] = [r0 * e for e in ratio]    # Annoying pythonic scalar * vector multiplication...
    
    with nazca.Cell(name=name) as mark:
        with nazca.Cell(instantiate=False) as pads:
            layerPolygon(layers=xs_pad, grow=grow_pad, poly= 
                             nazca.geom.circle(radius=r0, N=64)).put(0)
            layerPolygon(layers=xs_pad, grow=grow_pad, poly= 
                            nazca.geom.circle(radius=r1p, N=64) + \
                            nazca.geom.circle(radius=r1, N=64)[::-1]
                         ).put(0)
            layerPolygon(layers=xs_pad, grow=grow_pad, poly= 
                            nazca.geom.circle(radius=r2p, N=64) + \
                            nazca.geom.circle(radius=r2, N=64)[::-1]
                         ).put(0)
        pads.put(0)
        
        if (xs_background is not None):
            layerPolygon(layers=xs_background, poly=pads.bbox_polygon.tolist(), grow=grow_background).put(0)
        
        nazca.netlist.Annotation(layer='Annotation', text=( 'r0=%s, ratio=%s' % (str(r0), str(ratio)) )).put(-1.5*r2p, 1.1*r2p)
    
    return mark


# Circular TLM array [Pan 2013 https://doi.org/10/ghjtgf ]
def CircularTLM(xs_pad, xs_background=None, grow_pad=None, grow_background=10, r0=[5, 6.5, 9], ratio=[1, 1.5, 2.5], buffer=40, ratio_outer=1.25):
    """
    CircularTLM:
        Array of circular TLM patterns as in Pan 2013 ( https://doi.org/10/ghjtgf )

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
    """
    # Name and deduplicate
    name = 'TLM.Circular.'+str(r0)+'.'+str(ratio)+'.'+str(layerlist(xs_pad))
    name += '/'+str(layerlist(xs_background)) if xs_background is not None else ''
    if name in nazca.cfg.cellnames.keys():
        return nazca.cfg.cellnames[name]
    
    # Calculate spacing
    dxy = buffer + max(r0) * max(ratio) * ratio_outer * 1.5
    dx = (len(ratio)-1)*dxy + max(r0) * max(ratio) * ratio_outer + buffer*2
    dy = (len(r0)-1)*dxy + max(r0) * max(ratio) * ratio_outer + buffer*2
    
    polyOuter = [(-dx/2,dy/2), (dx/2,dy/2), (dx/2,-dy/2), (-dx/2,-dy/2)]
    polyBackground = polyOuter.copy()
    
    with nazca.Cell(name=name) as mark:
        with nazca.Cell(instantiate=False) as pads:
            y = (len(r0)-1)*dxy/2
            for r01 in r0:
                rOuter = r01 * max(ratio) * ratio_outer
                
                polyOuter += [(-dx/2, y + rOuter)]
                
                x = -(len(ratio)-1)*dxy/2
                for rat in ratio:
                    # Inner circle
                    layerPolygon(layers=xs_pad, grow=grow_pad, poly= 
                                    nazca.geom.circle(radius=r01*rat, N=64)
                                 ).put(x, y)
                    
                    # Construct outer pad boundary
                    polyOuter += addPolyOffset(reversed(nazca.geom.circle(radius=rOuter, N=64)), (x,y))
                    
                    x += dxy
                polyOuter += [(-dx/2, y + rOuter)]
                y -= dxy
        
        # Place inner pads
        pads.put(0, 0)
        
        # Place outer pad boundary
        layerPolygon(layers=xs_pad, grow=grow_pad, poly=polyOuter+[(-dx/2,dy/2)]).put(0)
        
        # Background, if any
        if (xs_background is not None):
            layerPolygon(layers=xs_background, poly=polyBackground, grow=grow_background).put(0)
        
        nazca.netlist.Annotation(layer='Annotation', text=( 'r0=%s, ratio=%s, ratio_outer=%s' % (str(r0), str(ratio), str(ratio_outer)) )).put(-dx/2, dy/2)
    
    return mark


# Layer label[s]
def LayerLabels(layers=None, size=25, labels=[]):
    """
    LayerLabels:
        Generates one or more labels for passed layers, in the form of '[#] layername'

    Args:
        layers (int, str, or list, optional): Layer(s) to generate labels for; if None, will use all layers
        size (int, optional): Text size (height). Defaults to 25.
        labels (list of str, optional): Labels to use instead of nazca-defined layer name

    Returns:
        Cell: Left-centered-origin cell of layer labels
    """
    # Process inputs
    if layers is None:
        layers = sorted(nazca.cfg.layer_table.layer)
    layers = flatten(layerlist(layers))
    labels = flatten(labels)
    
    # Name and deduplicate
    name = 'LayerLabels.'+str(layers)
    if name in nazca.cfg.cellnames.keys():
        return nazca.cfg.cellnames[name]
    
    # Iterate layers to get names and numbers
    for i in range(len(layers)):
        layers[i] = {'layer': layers[i], 'name': nazca.get_layer(layers[i])}
        if len(labels) > i:
            layers[i]['name'] = labels[i]
    
    # Build text objects
    dl = size/4
    dy = (len(layers)*size + (len(layers)-1)*dl)/2
    with nazca.Cell(name=name) as mark:
        for i in range(len(layers)):
            lLabel = '[%d] %s' % (layers[i]['layer'], layers[i]['name'])
            # Unless the text is rebuilt, bounding boxes are not computed!
            nazca.Font(fontfile='nazca').text(text=lLabel, height=size, 
                                              layer=layers[i]['layer'], align='lt').\
                                                  rebuild().put(0, dy - i*(size+dl))
    
    return mark


# Resolution test
## Single cell
def ResolutionTest(xs_pattern, xs_background=None, scale=1, 
                   grow_pattern=None, grow_background=5,
                   no_text=False):
    """ResolutionTest:
        Resolution test pattern including lines and checkerboard
    
    Args:
        xs_pattern (xs or layer): Pads defined in this cross-section [str] or layer(s) [int or list(int)].
        xs_background (xs or layer, optional): Background defined in this cross-section or layer.
        scale (float): Resolution to test in µm. Default 1.
        grow_* (int or list(int), optional): Grow associated geometry accordingly.
        no_text (bool): Omit text. Default False.
    
    Returns:
        Cell: Centered-origin resolution test pattern
    """
    boxSize = 20
    
    # Process inputs, name, and deduplicate
    xs_pattern = layerlist(xs_pattern)
    name = 'ResolutionTest.'+str(scale)+'.'+str(xs_pattern)
    if xs_background is not None:
        xs_background = layerlist(xs_background)
        name += '/'+str(xs_background)
    if name in nazca.cfg.cellnames.keys():
        return nazca.cfg.cellnames[name]
    
    # Define and deduplicate components to reduce filesize
    # Checkerboard component; smaller to define entire pattern
    checkerName = 'checkerboard.'+str([xs_pattern,grow_pattern])
    if checkerName in nazca.cfg.cellnames.keys():
        checkerboard = nazca.cfg.cellnames[checkerName]
    else:
        with nazca.Cell(instantiate=False) as checkerboard:
            box = layerRectangle(layers=xs_pattern, grow_layers=grow_pattern, dx=1)
            for j in range(7):
                for i in range(7):
                    if mod(i+j,2)==0:
                        box.put((i-7/2), (j-7/2))
        checkerboard = checkerboard.flatten(name=checkerName, instantiate=True)
    
    # Dot grid component; dot is a separate pattern as it has 24 vertices
    dotName = 'dot.'+str(xs_pattern)+'+'+str(grow_pattern)
    if dotName in nazca.cfg.cellnames.keys():
        dot = nazca.cfg.cellnames[dotName]
    else:
        with nazca.Cell(name=dotName) as dot:
            layerPolygon(layers=xs_pattern, 
                         poly=nazca.geom.circle(radius=(2**0.5)/2, N=24),
                         grow=grow_pattern).put(0)
    with nazca.Cell(instantiate=False) as dotgrid:
        for j in range(7):
            for i in range(7):
                if mod(i+j,2)==0:
                    dot.put(scale*(i-7/2), scale*(j-7/2), scale=scale)
    
    # Add checkerboards and lines
    with nazca.Cell(name=name, instantiate=True) as mark:
        if xs_background is not None:
            layerRectangle(layers=xs_background, grow_layers=grow_background, 
                           dx=scale*boxSize).put(0,0)
        
        # Checkerboard pattern
        checkerboard.put(scale*(-boxSize/2+7/2), 
                         scale*(-boxSize/2+7/2),
                         scale=scale)
        
        # Stripes
        for i in range(7):
            if mod(i,2)==1:
                # Horizontal stripes
                layerPolygon(layers=xs_pattern, 
                             poly=nazca.geom.rectangle(length=scale*(boxSize-7-i), 
                                                             height=scale, 
                                                             position=2), 
                             grow=grow_pattern).put(scale*(-boxSize/2+7-0.5),
                                                    scale*(-boxSize/2+i))
                # Vertical stripes
                layerPolygon(layers=xs_pattern, 
                             poly=nazca.geom.rectangle(length=scale*(boxSize-7-i), 
                                                             height=scale, 
                                                             position=2), 
                             grow=grow_pattern).put(scale*(boxSize/2-1-i),
                                                    scale*(-boxSize/2+i-0.5), 90)
        
        # Dot grid
        dotgrid.put(scale*(boxSize/2-7/2), 
                    scale*(boxSize/2-7/2))
        
        # Label
        if not no_text:
            cellShift(nazca.Font(fontfile='nazca').text(text='%.1f'%scale, 
                                                        height=scale*boxSize/2.5, 
                                                        layer=xs_pattern,
                                                        align='rb'),
                    ['bottom', 'right']).put(scale*(boxSize/2-9-1), 
                                            -scale*(boxSize/2-9))
    
    return mark


## Nested resolution tests
def matryoshkaResolution(xs_pattern, xs_background=None, 
                         grow_pattern=None, grow_background=5, 
                         scale_max=5, scale_min=0.2):
    """matryoshkaResolution:
        Nested, unlabeled resolution test patterns of reducing size
    
    Args:
        xs_pattern (xs or layer): Pads defined in this cross-section [str] or layer(s) [int or list(int)].
        xs_background (xs or layer, optional): Background defined in this cross-section or layer.
        grow_* (int or list(int), optional): Grow associated geometry accordingly.
        scale_max (float): Maximum test resolution [µm]
        scale_min (float): Minimum test resolution [µm]
    
    Returns:
        Cell: Centered-origin resolution test pattern
    """
    boxSize = 20    # from ResolutionTest
    
    # Process inputs, name, and deduplicate
    scale_max = round(scale_max, 1)
    xs_pattern = layerlist(xs_pattern)
    name = 'matryoshkaResolution.'+str(scale_max)+'-'+str(scale_min)+'.'+str(xs_pattern)
    if xs_background is not None:
        xs_background = layerlist(xs_background)
        name += '/'+str(xs_background)
    if name in nazca.cfg.cellnames.keys():
        return nazca.cfg.cellnames[name]
    
    # Generate mark
    with nazca.Cell(name=name, instantiate=True) as mark:
        if xs_background is not None:
            layerRectangle(layers=xs_background, 
                           dx=scale_max*boxSize+2*grow_background,
                           dy=scale_max*boxSize+2*grow_background, 
                           position=5).put(0,0)
        
        # Place test patterns
        x=0
        scale = []
        s = scale_max
        while s >= scale_min:
            scale += ['%.1f'%s]
            ResolutionTest(scale=s, xs_pattern=xs_pattern, xs_background=None,
                           grow_pattern=grow_pattern, no_text=True).put(-x, x)
            
            # Update variables
            x += s*boxSize/5
            s = round(s/2, 1)
        
        # Label
        nazca.Font(fontfile='nazca').text(text=', '.join(flip(scale)), strokewidth=(scale_max + 4*scale_min)/5, 
                                          layer=xs_pattern, align='cb').\
                                              put(0, scale_max*boxSize/2+5)
    
    return mark


## Improved test block
def ResolutionBlock(xs_pattern, xs_background=None, res=0.5, box_size=15, 
                    grow_pattern=None, grow_background=5):
    """ResolutionBlock:
        Resolution test pattern; all quadrants visible for test resolution, 
            diagonals for lower resolution, one side for uneven thresholding
    
    Args:
        xs_pattern (xs or layer): Cross-section [str] or layer(s) [int or list(int)] for test pattern.
        xs_background (xs or layer, optional): Cross-section or layer for background box.
        res (float): Resolution to test in µm.
        box_size (float): Test quadrant size in µm.
        grow_* (int or list(int), optional): Grow associated geometry accordingly.
    
    Returns:
        Cell: Centered-origin resolution test pattern
    """
    box_border = 0.5
    
    # Process inputs, name, and deduplicate
    xs_pattern = layerlist(xs_pattern)
    name = 'ResolutionBlock.'+str(xs_pattern)+'.'+str([res,box_size,grow_pattern])
    if xs_background is not None:
        xs_background = layerlist(xs_background)
        name += '/'+str(xs_background)+'.'+str(grow_background)
    if name in nazca.cfg.cellnames.keys():
        return nazca.cfg.cellnames[name]
    
    # Define and deduplicate basic pattern component to reduce filesize
    dotName = 'resDot.'+str(xs_pattern)+'+'+str(grow_pattern)
    if dotName in nazca.cfg.cellnames.keys():
        dot = nazca.cfg.cellnames[dotName]
    else:
        with nazca.Cell(name=dotName) as dot:
            layerPolygon(layers=xs_pattern, 
                         poly=nazca.geom.circle(radius=0.5, N=18),
                         grow=grow_pattern).put(0)
    
    # Define two basic patterns
    
    # # Checkerboard: larger circles for better thresholding balance determination
    # with nazca.Cell(instantiate=False) as checkerboard:
    #     # layerRectangle(layers=1000, dx=box_size, position=1).put(0,0)
    #     # layerRectangle(layers=1002, dx=box_size-2*box_border, position=1).put(box_border,box_border)
    #     dx = (box_size - 2*scale*(floor(Npt/2)+-0.5))/2
    #     dot.put(dx, dx, scale=scale*2**0.5, 
    #             array=[floor(Npt/2), [2*scale, 0], floor(Npt/2), [0, 2*scale]])
    #     dot.put(dx+scale, dx+scale, scale=scale*2**0.5, 
    #             array=[floor(Npt/2), [2*scale, 0], floor(Npt/2), [0, 2*scale]])
    
    # Grid: only resolves when resolution is at marked effective MFD
    with nazca.Cell(instantiate=False) as grid:
        scale = res*1.25   # Component scale needs to be 25% above resolution according to simulations
        Npt = max(1, floor((box_size-2*box_border)/scale))
        # layerRectangle(layers=1000, dx=box_size, position=1).put(0,0)
        # layerRectangle(layers=1002, dx=box_size-2*box_border, position=1).put(box_border,box_border)
        dot.put(-scale*(Npt-1)/2, -scale*(Npt-1)/2, scale=scale, 
                array=[Npt, [scale, 0], Npt, [0, scale]])
    
    # Lines: to determine sub-resolution patterning
    #   Place explicitly rather than arrayed for simplicity and avoiding instantiation
    with nazca.Cell(instantiate=False) as lines:
        # layerPolygon(layers=xs_pattern, poly=nazca.geom.frame(box_border, 
        #                                                       box_size-box_border, 
        #                                                       box_size-box_border)).\
        #                                                           put(-(box_size-2*box_border)/2, 
        #                                                               -(box_size-2*box_border)/2)
        Npt = max(1, floor(box_size/2/res))
        line = layerPolygon(layers=xs_pattern, grow=grow_pattern, 
                            poly=nazca.geom.rectangle(length=box_size, height=res, position=5))
        for n in range(Npt):
            line.put(0, -res*(Npt-1)+2*res*n)
    
    # Actual mark
    with nazca.Cell(name=name, instantiate=True) as mark:
        if xs_background is not None:
            layerRectangle(layers=xs_background, grow_layers=grow_background, 
                           dx=2*box_size+2*box_border).put(0,0)
        
        # Quadrants
        #   Not using flip/flop because arrays work incorrectly
        grid.put(-box_size/2, box_size/2)
        lines.put(box_size/2, box_size/2)
        grid.invert(bounds=nazca.geom.rectangle(length=box_size, height=box_size, position=5),
                    mod_layers=xs_pattern, instantiate=False).put(-box_size/2, -box_size/2)
        lines.put(box_size/2, -box_size/2, 90)
        # checkerboard.invert(bounds=nazca.geom.rectangle(length=box_size, height=box_size),
                            # mod_layers=xs_pattern, instantiate=False).put(0, -box_size)
        
        # Border for quick visual identification
        layerPolygon(layers=xs_pattern, 
                     poly=nazca.geom.frame(2*box_border, 2*box_size+2*box_border, 2*box_size+2*box_border)).\
            put(-box_size, -box_size)
    
    return mark


# Set of multiple resolution blocks, with labels
def ResolutionBlockSet(xs_pattern, res=[0.6, 0.7, 0.8, 0.9, 1.0, 1.1], 
                       box_size=15, with_label=True, label_size=20, 
                       xs_background=None, grow_pattern=None, grow_background=5):
    """ResolutionBlockSet:
        Set of resolution block test patterns at different resolutions, with labels, in two rows
    
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
    """
    # Process inputs, name, and deduplicate
    xs_pattern = layerlist(xs_pattern)
    name = 'ResolutionBlockSet.'+str(xs_pattern)+'.'+str(res)+'.'+str([box_size,with_label,label_size,grow_pattern])
    if xs_background is not None:
        xs_background = layerlist(xs_background)
        name += '/'+str(xs_background)+'.'+str(grow_background)
    if name in nazca.cfg.cellnames.keys():
        return nazca.cfg.cellnames[name]
    
    # Resolution test in two rows
    with nazca.Cell(name=name, instantiate=True) as mark:
        dx = max(box_size*4/3, label_size*2) + 4
        x = -dx*(len(res)/2-1)/2
        y = dx/2
        for r in res:
            # print(str([x,y]))
            ResolutionBlock(xs_pattern=xs_pattern, xs_background=xs_background, res=r, 
                            box_size=box_size, grow_pattern=grow_pattern, 
                            grow_background=grow_background).put(x, y)
            if with_label:
                label(text='%.2g' % r, height=label_size, layer=xs_pattern, grow=grow_pattern, 
                      origin=['center', 'lower' if y >= 0 else 'upper']).\
                          put(x, y + (dx/2 if y >= 0 else -dx/2))
            x += (dx if y<0 else 0)
            y = mod(y+1.5*dx, 2*dx) - dx/2
    return mark
