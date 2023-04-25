Module marks_NanoFab.stepper2
=============================
Stepper2 (GCA AutoStep 200) tool markings.
(c) Michael Nickerson 2023

Functions
---------

    
`GlobalAlign(layer=1, background=None, grow_background=10)`
:   Global alignment mark for AutoStep 200 ('star').
    
    Args:
        layer (xs or layer): Cross-section [str] or layer(s) to place mark in
        background (xs or layer, optional): Bounding background defined in this cross-section or layer
        grow_background (int or list(int), optional): Grow associated geometry accordingly
    
    Returns:
        Cell: Centered-origin global alignment mark

    
`LocalAlign(layer=1, background=None, grow_background=10)`
:   Local DFAS alignment mark for AutoStep 200.
    
    Args:
        layer (xs or layer): Cross-section [str] or layer(s) to place mark in
        background (xs or layer, optional): Bounding background defined in this cross-section or layer
        grow_background (int or list(int), optional): Grow associated geometry accordingly
    
    Returns:
        Cell: Centered-origin local alignment mark

    
`ReticleMarks(layer=1, scale=1, darkfield=True, inner=[-9200, -9200, 9200, 9200], innerstepback=0)`
:   Reticle alignment marks for AutoStep 200.
    
    Args:
        layer (xs or layer): Cross-section [str] or layer(s) to place reticle
        scale (1 or 5): mask (1) or reticle (5) scale
        darkfield (bool): mask polarity; if true, geometry = clear
        inner (list(float)): lower-left and upper-right coordinates of inner boundary for lightfield inversion
        innerstepback (float): stepback from inner boundary for lightfield inversion
    
    Returns:
        Cell: Centered-origin reticle markings

    
`addReticles(cell, layers=None, darkfield=True, instantiate=True)`
:   Add 5-inch reticle(s) for AutoStep 200 to a cell.
    
    Args:
        cell (Cell): cell to place in middle of mask at (0,0)
        layers (int or str): Layers to build reticles for; default all present
        darkfield (bool or list(bool)): mask polarity; if true, geometry = clear; list corresponds to <layers>
    
    Returns:
        Cell: New cell with proper reticles in MASK COORDINATES

    
`die(layer=1001, size=[14800, 14800])`
:   Generates die outline for AutoStep 200.
    
    Args:
        layer (int, optional): Outline layer; likely 'comment' layer. Defaults to 1001.
        size (list, optional): Die size, maximum [14800, 14800]. Defaults to maximum.
    
    Returns:
        Cell: Centered-origin outline