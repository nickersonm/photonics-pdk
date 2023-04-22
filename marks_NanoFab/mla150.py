#!/usr/bin/env python3
# 
# Provides standard UCSB NanoFab tool markings
#   This file for Heidelberg MLA150
#       https://wiki.nanotech.ucsb.edu/wiki/Maskless_Aligner_(Heidelberg_MLA150)
# 
# @Authors: Michael Nickerson
# @email: nickersonm@ece.ucsb.edu
# 2021(c)
# 

"""
Heidelberg MLA150 tool markings
"""

import nazca
from .fCommon import layerlist, layerPolygon, flatten

# High-mag alignment mark
#   https://wiki.nanotech.ucsb.edu/wiki/MLA150_-_CAD_Files_and_Templates
def AlignHighMag(layer=1, background=None, grow_background=10):
    """
    AlignHighMag High-magnification alignment mark for Heidelberg MLA150 ('cross')

    Args:
        layer (xs or layer): Cross-section [str] or layer(s) to place mark in
        background (xs or layer, optional): Bounding background defined in this cross-section or layer
        grow_background (int or list(int), optional): Grow associated geometry accordingly

    Returns:
        Cell: Centered-origin alignment mark
    """
    lineWidth=8
    size=150
    
    layer = flatten(layerlist(layer))
    
    markName = 'MLA150.AlignHighMag.'+str(layer)
    markName += '/'+str(layerlist(background)) if background is not None else ''
    if markName in nazca.cfg.cellnames.keys():
        return nazca.cfg.cellnames[markName]
    
    with nazca.Cell(instantiate=False) as mark:
        # layerPolygon(layers=layer, 
        #              poly=[nazca.geom.rectangle(length=size, height=line, position=5),
        #                    nazca.geom.rectangle(length=line, height=size, position=5)]).put(0)
        with nazca.Cell(instantiate=False) as arm:
            layerPolygon(layers=layer,
                         poly=[[0,0], [lineWidth, lineWidth/2], [size/2, lineWidth/2],
                               [size/2, -lineWidth/2], [lineWidth, -lineWidth/2]]).put(0)
        for a in [0, 90, 180, 270]:
            arm.put(0,0,a)
        
        if (background is not None):
            layerPolygon(layers=background, 
                         poly=nazca.geom.rectangle(length=size, height=size, position=5), 
                         grow=grow_background).put(0)
    
    return mark.flatten(name=markName, instantiate=True)

