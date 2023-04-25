Module marks_NanoFab.fCommon
============================
Common utility functions.
(c) Michael Nickerson 2023

Functions
---------

    
`addPolyOffset(poly, offset)`
:   Add offset to all points in a polygon.

    
`cellFlipPins(cell, instantiate=None)`
:   Flip cell pins - only applies to a0, b0.

    
`cellHeight(cell)`
:   Find the height of a cell (y-span).

    
`cellPolygons(cell, layers)`
:   Return all polygons in <layers> as point lists.
    
    Args:
        cell (Cell): cell to get polygons from
        layers (int or str): names of layers to retreive polygons from
    
    Returns:
        List: list of polygon points

    
`cellShift(cell, loc='center')`
:   Adjust cell origin to bounding box center or given edge or corner.

    
`cellWidth(cell)`
:   Find the width of a cell (x-span).

    
`flat(l)`
:   

    
`flatten(l)`
:   Flatten a list.

    
`layerPolygon(layers, poly, grow=None, jointype='miter')`
:   Build a cell out of a polygon on a given set of layers (or cross-section).

    
`layerRectangle(layers, dx, dy=None, offset=(0, 0), grow_layers=None, position=5)`
:   Build a rectangle in a given cross-section.

    
`layerlist(xs_or_layer)`
:   Return a list of layers when passed an xs [cross-section] or layer number or list of layer numbers.

    
`putCorners(cell, diesize, inset=[50, 50], shift=True, quadrants=[1, 2, 3, 4], flip=False, flop=False)`
:   Place given cell in corners of given diesize, distributed along x by 'spacing'.
    
    Args:
        cell (Cell): cell to place
        diesize (list<float>): size [dx,dy] for centered die, or lower-left and upper-right corner [[x0, y0], [x1, y1]]
        inset (float): xy or [x,y] spacing from edge
        shift (bool): shift origin of cell to cell corners before placing
        quadrants (list<float>): place in these quadrants
        flip, flop (bool): flip/flop cells at corners for symmetry; does not adjust cell origin
    
    Returns:
        List: list of polygon points