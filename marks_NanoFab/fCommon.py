#!/usr/bin/env python3
# 
# Provides common utility functions
# 
# @Authors: Michael Nickerson
# @email: nickersonm@ece.ucsb.edu
# 2020(c)
# 

"""
Common utility functions.
(c) Michael Nickerson 2023
"""

import nazca


# Return a list of layers when passed an xs [cross-section] or layer number or list of layer numbers
def layerlist(xs_or_layer):
    """Return a list of layers when passed an xs [cross-section] or layer number or list of layer numbers."""
    ll = []
    if isinstance(xs_or_layer, list):
        # List of cross-sections or layer numbers
        for l in xs_or_layer:
            lx = layerlist(l)
            if lx is not None:
                ll.append(lx if len(lx) > 1 else lx[0])
    elif isinstance(xs_or_layer, str) and (xs_or_layer in nazca.cfg.XSdict):
        # Is a valid cross-section
        for l in nazca.mask_elements.layeriter(xs=xs_or_layer):
            lx = layerlist(l[0])
            if lx is not None:
                ll.append(lx if len(lx) > 1 else lx[0])
    else:
        if xs_or_layer in nazca.mask_layers.get_layers()['layer'] or \
            xs_or_layer in nazca.mask_layers.get_layers()['layer'].values:
            # Attempt to validate
            ll.append(nazca.get_layer_tuple(layer=xs_or_layer).layer)
        elif isinstance(xs_or_layer, float) or isinstance(xs_or_layer, int):
            # Just add as is - will force a new layer
            ll.append(xs_or_layer)
    
    return list(flatten(ll))


# Build a cell out of a polygon on a given set of layers (or cross-section)
def layerPolygon(layers, poly, grow=None, jointype='miter'):
    """Build a cell out of a polygon on a given set of layers (or cross-section)."""
    layers = layerlist(layers)
    if not isinstance(grow, list):
        grow = len(layers) * ([0] if grow is None else [grow])
    elif isinstance(grow, list) and (len(grow) != len(layers)):
        nazca.main_logger(msg='Growth specified as list but not the same length as layers!', level='warning')
    if len(poly[0]) < 3:
        poly = [poly]
    
    with nazca.Cell(instantiate=False) as polyCell:
        for i in range(len(layers)):
            if grow[i] != 0:
                newPoly = [p for p in nazca.clipper.grow_polygons(poly, grow[i], accuracy=1e-4, jointype=jointype)]
                newPoly = [p for p in nazca.clipper.merge_polygons(newPoly, accuracy=1e-4)]
                # newPoly = [p.points for p in
                #            gdstk.offset(polygons=poly, distance=grow[i], use_union=True)]
            else:
                newPoly = poly
            
            if len(newPoly[0]) < 3:
                newPoly = [newPoly]
            for p in newPoly:
                nazca.Polygon(points=p, layer=layers[i]).put()
    return polyCell


# Build a rectangle in a given cross-section
def layerRectangle(layers, dx, dy=None, offset=(0,0), grow_layers=None, position=5):
    """Build a rectangle in a given cross-section."""
    dy = dy if dy is not None else dx
    return layerPolygon(layers=layers, poly=nazca.geometries.rectangle(length=dx, height=dy, position=position, shift=offset), grow=grow_layers)

# Flatten a list
def flat(l):
    for t in l:
        if hasattr(t, '__iter__') and not isinstance(t, str):
            yield from flat(t)
        else:
            yield t
def flatten(l):
    """Flatten a list."""
    return list(flat(l))


## Various simple calculations
# Pin separation
def dx(cell):
    """x-separation between input and output pins."""
    return round(abs(cell.pinout.x - cell.pinin.x), 3)
def dy(cell):
    """y-separation between input and output pins."""
    return round(abs(cell.pinout.y - cell.pinin.y), 3)

# Add offset to a polygon tuple
def addPolyOffset(poly, offset):
    """Add offset to all points in a polygon."""
    return [(p[0] + offset[0], p[1] + offset[1]) for p in poly]

# Width of a cell
def cellWidth(cell):
    """Find the width of a cell (x-span)."""
    if cell.bbox is None:   # May not have bbox computed
        cell = cell.rebuild(instantiate=False)
    return abs(cell.bbox[2] - cell.bbox[0])

# Height of a cell
def cellHeight(cell):
    """Find the height of a cell (y-span)."""
    if cell.bbox is None:   # May not have bbox computed
        cell = cell.rebuild(instantiate=False)
    return abs(cell.bbox[3] - cell.bbox[1])

# Change offset of cell
def cellShift(cell, loc='center'):
    """Adjust cell origin to bounding box center or given edge or corner."""
    if cell.bbox is None:   # May not have bbox computed
        cell = cell.rebuild(instantiate=False)
    
    os = [0,0]
    
    if 'center' in loc or 'middle' in loc:
        os = [ -(cell.bbox[0] + cell.bbox[2])/2,
               -(cell.bbox[1] + cell.bbox[3])/2 ]
    
    if sum([c in ['left', 0] for c in loc]) > 0:
        os[0] = -cell.bbox[0]
    
    if sum([c in ['right', 2] for c in loc]) > 0:
        os[0] = -cell.bbox[2]
    
    if sum([c in ['bottom', 'lower', 1] for c in loc]) > 0:
        os[1] = -cell.bbox[1]
    
    if sum([c in ['top', 'upper', 3] for c in loc]) > 0:
        os[1] = -cell.bbox[3]
    
    # More elegant except for handling 'center'
    # os =  [ -cell.bbox[0] * sum([c == 'left' or c == 0 for c in loc]) +
    #         -cell.bbox[2] * sum([c == 'right' or c == 2 for c in loc]),
    #         -cell.bbox[1] * sum([c == 'bottom' or c == 1 for c in loc]) +
    #         -cell.bbox[3] * sum([c == 'top' or c == 3 for c in loc])
    #        ]
    
    with nazca.Cell(instantiate=False) as C:
        cell.put(os[0], os[1])
    
    return C

# Flip cell pins - only applies to a0, b0
def cellFlipPins(cell, instantiate=None):
    """Flip cell pins - only applies to a0, b0."""
    name = cell.cell_name + '_flip'
    instantiate = cell.instantiate
    if name in nazca.cfg.cellnames.keys() and instantiate == True:
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(instantiate=instantiate, name=name) as C:
        pinsin = list(cell.pin.keys())
        pinsout = pinsin.copy()
        ia, ib = pinsin.index('a0'), pinsin.index('b0')
        pinsout[ia], pinsout[ib] = 'b0', 'a0'
        cell.put(0).raise_pins(namesin=pinsin, namesout=pinsout)
        C.length_geo = cell.length_geo if hasattr(cell, 'length_geo') else (dx(cell)**2 + dy(cell)**2)**0.5
    return C

# Place cells in corners
def putCorners(cell, diesize, inset=[50,50], 
               shift=True, quadrants=[1,2,3,4],
               flip=False, flop=False):
    """Place given cell in corners of given diesize, distributed along x by 'spacing'.
    
    Args:
        cell (Cell): cell to place
        diesize (list<float>): size [dx,dy] for centered die, or lower-left and upper-right corner [[x0, y0], [x1, y1]]
        inset (float): xy or [x,y] spacing from edge
        shift (bool): shift origin of cell to cell corners before placing
        quadrants (list<float>): place in these quadrants
        flip, flop (bool): flip/flop cells at corners for symmetry; does not adjust cell origin
    
    Returns:
        List: list of polygon points
    """
    # Process input
    if len(diesize) == 1:   # Legacy square input
        diesize = [diesize, diesize]
    if isinstance(diesize[0], int):
        diesize = [[-diesize[0]/2, -diesize[1]/2], [diesize[0]/2, diesize[1]/2]]
    if not isinstance(inset, list):
        inset = [inset, inset]
    if shift:
        os = ['left', 'right', 'bottom', 'top']
    else:
        os = [None, None, None, None]
    
    if 3 in quadrants:  # Lower-left
        cellShift(cell, [os[0], os[2]]).\
            put(diesize[0][0] + inset[0], diesize[0][1] + inset[1])
    if 2 in quadrants:  # Upper-left
        cellShift(cell, [os[0], os[3-int(flip)]]).\
            put(diesize[0][0] + inset[0], diesize[1][1] - inset[1], flip=flip)
    if 1 in quadrants:  # Upper-right
        cellShift(cell, [os[1-int(flop)], os[3-int(flip)]]).\
            put(diesize[1][0] - inset[0], diesize[1][1] - inset[1], flop=flop, flip=flip)
    if 4 in quadrants:  # Lower-right
        cellShift(cell, [os[1-int(flop)], os[2]]).\
            put(diesize[1][0] - inset[0], diesize[0][1] + inset[1], flop=flop)


# Get polygons of a cell and layer
def cellPolygons(cell, layers):
    """Return all polygons in <layers> as point lists.
    
    Args:
        cell (Cell): cell to get polygons from
        layers (int or str): names of layers to retreive polygons from
    
    Returns:
        List: list of polygon points
    """
    layers = layerlist(layers)
    polys = []
    for params in nazca.cell_iter(cell, flat=True):
        if params.cell_start:
            for poly, xy, _ in params.iters['polygon']:
                if len(set(layerlist(poly.layer)) & set(layers)) > 0:
                    polys.append(xy)
    return polys
    