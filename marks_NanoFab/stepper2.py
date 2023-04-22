#!/usr/bin/env python3
# 
# Provides standard UCSB NanoFab tool markings
#   This file for GCA AutoStep 200
# 
# @Authors: Michael Nickerson
# @email: nickersonm@ece.ucsb.edu
# 2020(c)
# 

"""
Stepper2 (GCA AutoStep 200) tool markings
"""

import nazca
from . import pathGDS
from .fCommon import layerlist, layerRectangle, layerPolygon, flatten, cellPolygons, cellShift
from .utility import die as utility_die

from marks_NanoFab import fCommon


# Global alignment mark
def GlobalAlign(layer=1, background=None, grow_background=10):
    """
    GlobalAlign Global alignment mark for AutoStep 200 ('star')

    Args:
        layer (xs or layer): Cross-section [str] or layer(s) to place mark in
        background (xs or layer, optional): Bounding background defined in this cross-section or layer
        grow_background (int or list(int), optional): Grow associated geometry accordingly

    Returns:
        Cell: Centered-origin global alignment mark
    """    
    layer = flatten(layerlist(layer))
    
    markName = 'S2.GlobalAlign.'+str(layer)
    markName += '/'+str(layerlist(background)) if background is not None else ''
    if markName in nazca.cfg.cellnames.keys():
        return nazca.cfg.cellnames[markName]
    
    S2_GlobalAlign = nazca.cfg.cellnames['S2.GlobalAlign'] if 'S2.GlobalAlign' in nazca.cfg.cellnames.keys() else \
                     nazca.load_gds(filename=pathGDS+'GlobalAlign.gds', instantiate=False, 
                                    cellname='GlobalAlign', newcellname='S2.GlobalAlign', scale=1)
    
    with nazca.Cell(name=markName) as mark:
        for l in layer:
            S2_GlobalAlign.rebuild(flat=True, layermapmode='none', 
                                   layermap={1:l}, instantiate=False).put(0)
        if (background is not None):
            layerPolygon(layers=background, poly=S2_GlobalAlign.bbox_polygon.tolist(), grow=grow_background).put(0)
    
    return mark


# Local alignment mark
def LocalAlign(layer=1, background=None, grow_background=10):
    """
    LocalAlign Local DFAS alignment mark for AutoStep 200
    
    Args:
        layer (xs or layer): Cross-section [str] or layer(s) to place mark in
        background (xs or layer, optional): Bounding background defined in this cross-section or layer
        grow_background (int or list(int), optional): Grow associated geometry accordingly

    Returns:
        Cell: Centered-origin local alignment mark
    """    
    layer = flatten(layerlist(layer))
    
    markName = 'S2.LocalAlign.'+str(layer)
    markName += '/'+str(layerlist(background)) if background is not None else ''
    if markName in nazca.cfg.cellnames.keys():
        return nazca.cfg.cellnames[markName]
    
    S2_LocalAlign = nazca.cfg.cellnames['S2.LocalAlign'] if 'S2.LocalAlign' in nazca.cfg.cellnames.keys() else \
                    nazca.load_gds(filename=pathGDS+'LocalAlign.gds', instantiate=False, 
                                   cellname='LocalAlign', newcellname='S2.LocalAlign', scale=1)
    
    with nazca.Cell(name=markName) as mark:
        for l in layer:
            S2_LocalAlign.rebuild(flat=True, layermapmode='none', 
                                  layermap={1:l}, instantiate=False).put(0)
        if (background is not None):
            layerPolygon(layers=background, poly=S2_LocalAlign.bbox_polygon.tolist(), grow=grow_background).put(0)
    
    return mark


# Die definition
def die(layer=1001, size=[14800, 14800]):
    """
    die Generates die outline for AutoStep 200

    Args:
        layer (int, optional): Outline layer; likely 'comment' layer. Defaults to 1001.
        size (list, optional): Die size, maximum [14800, 14800]. Defaults to maximum.
    
    Returns:
        Cell: Centered-origin outline
    """
    size = [min(14800, size[0]), min(14800, size[1])]
    # Check size
    if size[0] > 14800 or size[1] > 14800:
        nazca.main_logger(
            msg='stepper2.die(): die exceeds 14.8x14.8mm!  May not expose properly.', 
            level='warning')
    
    return utility_die(size, layer)


# Mask marking
def ReticleMarks(layer=1, scale=1, darkfield=True, 
                 inner=[-9200, -9200, 9200, 9200], innerstepback=0):
    """
    Reticle alignment marks for AutoStep 200
    
    Args:
        layer (xs or layer): Cross-section [str] or layer(s) to place reticle
        scale (1 or 5): mask (1) or reticle (5) scale
        darkfield (bool): mask polarity; if true, geometry = clear
        inner (list(float)): lower-left and upper-right coordinates of inner boundary for lightfield inversion
        innerstepback (float): stepback from inner boundary for lightfield inversion
    
    Returns:
        Cell: Centered-origin reticle markings
    """
    layer = flatten(layerlist(layer))
    
    # Deduplicate
    markName = 'S2.ReticleMark.'+str(layer)+'.'+ ('dark' if darkfield else 'light')
    if markName in nazca.cfg.cellnames.keys():
        return nazca.cfg.cellnames[markName]
    
    # Load
    S2_Reticle = nazca.cfg.cellnames['S2.ReticleAlign'] if 'S2.ReticleAlign' in nazca.cfg.cellnames.keys() else \
                 nazca.load_gds(filename=pathGDS+'Stepper2_Reticle_Dark.gds', instantiate=False, 
                                cellname='Reticle_Align', newcellname='S2.LocalAlign', scale=1)
    reticlePolys = cellPolygons(S2_Reticle, 1)
    
    if not darkfield:
        # Need to invert the mask for the outer boundary
        darkPoly = nazca.clipper.merge_polygons(
            [nazca.geom.rectangle(length=12700*2, height=12700*2, position=5),
             [[inner[0] - innerstepback, inner[1] - innerstepback],
              [inner[0] - innerstepback, inner[3] + innerstepback],
              [inner[2] + innerstepback, inner[3] + innerstepback],
              [inner[2] + innerstepback, inner[1] - innerstepback]] ], accuracy=1e-4)
        reticlePolys = nazca.clipper.diff_polygons(darkPoly, reticlePolys, accuracy=1e-4)
    
    with nazca.Cell(name=markName) as mark:
        for l in layer:
            for p in reticlePolys:
                nazca.Polygon(points=p, layer=l).put(0, scale=scale)
    
    return mark


# Prepare a reticle with marks and cell
def addReticles(cell, layers=None, darkfield=True, instantiate=True):
    """
    Add 5-inch reticle(s) for AutoStep 200 to a cell
    
    Args:
        cell (Cell): cell to place in middle of mask at (0,0)
        layers (int or str): Layers to build reticles for; default all present
        darkfield (bool or list(bool)): mask polarity; if true, geometry = clear; list corresponds to <layers>
    
    Returns:
        Cell: New cell with proper reticles in MASK COORDINATES
    """
    # TODO: quad mask?
    
    # Check cell size
    if (cell.bbox[2] - cell.bbox[0] > 14800) or (cell.bbox[3] - cell.bbox[1] > 14800):
        nazca.main_logger(
            msg='writeReticles: cell "%s" exceeds 14.8x14.8mm!  May not expose properly.' % (cell.cell_name), 
            level='warning')
    
    # Determine layers
    if layers is None:
        layers = []
        for params in nazca.cell_iter(cell, flat=True):
            for poly, _, _ in params.iters['polygon']:
                layers.append(poly.layer)
        layers = list(set(layers))   # Reduce to unique values
    else:
        layers = flatten(layerlist(layers))
    
    # Resize darkfield if needed
    if not isinstance(darkfield, list):
        darkfield = [darkfield]
    if len(darkfield) != len(layers):
        for _ in range(len(layers)-len(darkfield)):
            darkfield.append(darkfield[-1])
    
    # Generate name
    cellName = cell.cell_name + '.Reticles'
    
    # Get reticle alignment marks
    marks = []
    for layer, dark in zip(layers, darkfield):
        align = ReticleMarks(layer=layer, darkfield=dark, inner=cell.bbox, innerstepback=20)
        marks.append(align)
    
    # Assemble reticle markings
    with nazca.Cell(name=cellName, instantiate=instantiate) as reticle:
        # Place reticle alignment marks
        for align in marks:
            align.put(0)
        
        # Place cell with given layers
        cell.put(0)
    
    return reticle
