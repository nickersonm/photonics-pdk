#!/usr/bin/env python3
# 
# Provides utility functions for pdk_Fab6
# 
# @Authors: Michael Nickerson
# @email: nickersonm@ece.ucsb.edu
# 2023(c)
# 

"""
Utility functions for pdk_Fab6.
(c) Michael Nickerson 2023
"""

import nazca

from marks_NanoFab.utility import layerlist
from marks_NanoFab.fCommon import cellHeight, cellWidth, cellShift, cellFlipPins

from math import sin, pi

# For better polygon handling in some places
import gdstk    # https://heitzmann.github.io/gdstk/reference_python.html



### General helper functions
# Typical safe taper length for deep etched ridges from simulation
taperLen = lambda w1, w2 : round(abs(w2-w1)*30, 1)
"""Typical safe taper length for deep etched ridges from simulation"""

# Flip a cell
nazca.Cell.flip = lambda self : cellFlipPins(self)
nazca.Cell.reverse = nazca.Cell.flip
"""Flip a cell with `nazca.cell.flip`."""

# Change a cell name
def cellChangeName(cell, newname):
    """Change a cell name."""
    nazca.cfg.cellnames.pop(cell.cell_name)
    cell.cell_name = newname
    nazca.cfg.cellnames[newname] = cell

# Parse a pin/location/node into a real Pin object, with all transformations applied
def parsePin(pin, default='out', rot=False):
    """Parse a pin/location/node into a real Pin object, with all transformations applied."""
    pin, T = nazca.netlist.parse_pin(pin, default=default)
    if T != [0, 0, 0] or not isinstance(pin, nazca.Node):
        pin  = pin.move(*T)
        if default=='in':
            pin = pin.rot(180)
    return pin

# Get the geometric length between two pins with `nazca.cell.geolen`
def pinGeoDiff(p1, p2):
    """Get the geometric length between two pins with `nazca.cell.geolen`."""
    return (nazca.diff(p1, p2)[0]**2 +nazca.diff(p1, p2)[1])**0.5
nazca.Cell.geolen = lambda self : pinGeoDiff(self.pinin, self.pinout)

# Simple debugging function to plot polygons
def plotPolygon(p):
    """Simple debugging function to plot polygons."""
    with nazca.Cell(instantiate=False) as tmp:
        for p in p:
            nazca.Polygon(points=p, layer=1001).put(0)
    nazca.export_plt(tmp, clear=False)
    return None

# Generic array of cells
def cellArray(cells, space=20, nx=3):
    """Generic array of cells."""
    nx = round(nx)
    # Break into chunks; alternate: numpy.array_split
    cells = [cells[i:i + nx] for i in range(0, len(cells), nx)]
    
    # Place cells
    x = 0
    y = 0
    with nazca.Cell(instantiate=False) as a:
        for g in cells:
            with nazca.Cell(instantiate=False) as r:
                for c in g:
                    c.put(x, y)
                    x += cellWidth(c) + space
            cellShift(r).put(0, y, 'org')
            x = 0
            y -= cellHeight(r) + space
    
    return cellShift(a)

# Get polygons from a cell layer
def cell_get_polygons(self, layers):
    """Return all polygons in <layers> as point lists.
    Also callable as `nazca.cell.get_polygons`.
    
    Args:
        layers (list of str): names of layers to retreive polygons from
    
    Returns:
        List: list of polygon points
    """
    layers = layerlist(layers)
    polys = []
    for params in nazca.cell_iter(self, flat=True):
        if params.cell_start:
            for poly, xy, _ in params.iters['polygon']:
                if len(set(layerlist(poly.layer)) & set(layers)) > 0:
                    polys.append(xy)
                    # TODO: this might not get relative shifts in all cases?
    
    return polys
nazca.Cell.get_polygons = cell_get_polygons

# Remove a layer from a cell
def cell_remove_layer(self, layers, name=None, instantiate=None):
    """Remove all polygons in `layers`.
    Also callable as `nazca.cell.remove_layer`.
    
    Args:
        layers (list of str or int or xs): names of layers to delete polygons from
        name (str): new name for the cell
        instantiate (bool): instantiate the new cell?
    
    Returns:
        Cell: new cell with removed polygons
    """
    layers = layerlist(layers)
    
    # Deduplicate
    if name is None:
        name = self.cell_name+'-'+str(layers)
        if name in nazca.cfg.cellnames.keys():
            return nazca.cfg.cellnames[name]
    
    if instantiate is None:
        instantiate = self.instantiate
    
    # Rebuild without desired layer
    layers = {l:None for l in layers}
    c = self.rebuild(flat=True, instantiate=False, 
                     layermap=layers, layermapmode='all', 
                     infolevel=0)
    
    # TODO: Copy pins; this probably trips up pin connections?
    c.pin = self.pin
    c.default_in = self.default_in
    c.default_out = self.default_out
    
    cellChangeName(c, name)
    c.instantiate = instantiate
    c.length_geo = self.length_geo if hasattr(self, 'length_geo') else None
    c.openlen = self.openlen if hasattr(self, 'openlen') else None
    
    return c
nazca.Cell.remove_layer = cell_remove_layer

# Flatten a cell, including merging all polygons per layer
#   Don't use when the cell includes holes or concentric items
def cell_flatten(self, name=None, instantiate=None):
    """Rebuild cell to merge all polygons per layer.
    Also callable as `nazca.cell.flatten`.
    
    Args:
        name (str): new name for the cell
        instantiate (bool): instantiate the new cell?
    
    Returns:
        Cell: new cell with merged polygons
    """
    
    # Deduplicate
    if name is None:
        name = self.cell_name+'_flat'
    if name in nazca.cfg.cellnames.keys() and nazca.cfg.cellnames[name].flattened:
        return nazca.cfg.cellnames[name]
    
    if instantiate is None:
        instantiate = self.instantiate
    
    # Get relevant layers and polygons
    layers = []
    polys = nazca.defaultdict(list)
    for params in nazca.cell_iter(self, hierarchy='flat'):
        for poly, xy, _ in params.iters['polygon']:
            polys[poly.layer].append(xy)
            layers.append(poly.layer)
    layers = list(set(layers))   # Reduce to unique values
    
    # Merge polygons per layer
    for lay, xy in polys.items():
        # Merge
        if len(xy) > 0:
            polys[lay] = [p.points for p in 
                          gdstk.boolean(xy, xy, 'or', precision=1e-3) ]
        
        # newpolys=[]
        # # Grow very slightly to avoid imperfect overlaps
        # for poly in xy:
        #     if not isinstance(poly, list):
        #         poly = poly.tolist()
        #     if poly:
        #         for g in nazca.clipper.\
        #             grow_polygons([poly], grow=1e-3, accuracy=1e-4, jointype='square'):
        #                 newpolys.append(g)
        
        # # Merge
        # if len(xy) > 0:
        #     polys[lay] = []
        #     for poly in nazca.clipper.merge_polygons(newpolys, accuracy=1e-4):
        #         polys[lay].append(poly)
    
    # New cell
    with nazca.Cell(name=name, instantiate=instantiate) as c:
        # Add polygons
        for lay, xy in polys.items():
            for poly in xy:
                nazca.Polygon(points=poly, layer=lay).put(0)
    
    # Copy pins
    c.pin = self.pin
    c.pinin = self.pinin
    c.pinout = self.pinout
    
    # Copy other items
    c.length_geo = self.length_geo if hasattr(self, 'length_geo') else None
    c.flattened = True
    
    return c
nazca.Cell.flatten = cell_flatten

# Diff/clip layers of a cell
def cell_polydiff_layer(self, clip_poly, mod_layer=1001, operation='not'):
    """Rebuild cell to clip all polygons in specified layer(s); implicitly flattens mod_layer.
    Also callable as `nazca.cell.polydiff_layer`.
    
    Args:
        clip_poly (list): polygon geometry to use as clip
        mod_layer (int or str): layer(s) to modify; mod_layer' = mod_layer - clip_poly
    
    Returns:
        Cell: existing cell with clipped polygons
    """
    layers = layerlist(mod_layer)
    
    # # Grow and merge clip_poly
    # clip_poly = [ nazca.clipper.
    #     grow_polygons([poly], grow=1e-3, accuracy=1e-4, jointype='miter')[0]\
    #         for poly in clip_poly ]
    # clip_poly = nazca.clipper.merge_polygons(clip_poly, accuracy=1e-4)
    
    # Make sure clip_poly is 1-depth list of points
    clip_poly = [clip_poly]   # Should be a list of polygons
    while isinstance(clip_poly[0][0][0], list) or isinstance(clip_poly[0][0][0], tuple):
        clip_poly = clip_poly[0]
    
    # Get relevant polygons in layer(s)
    polys = nazca.defaultdict(list)
    for lay in layers:
        laypoly = self.get_polygons(lay)
        if len(laypoly) > 0:
            polys[lay] = laypoly
    if len(polys) < 1:
        return self # Nothing to do
    
    for params in nazca.cell_iter(self, flat=True):
        if params.cell_start:
            keepPolys = []
            for poly in params.cell.polygons:
                if len(set(layerlist(poly[1].layer)) & set(layers)) == 0:
                    # Keep polygons not in mod_layer
                    keepPolys.append(poly)
            params.cell.polygons = keepPolys
    
    # Clip and add mod_layer polys
    for lay, xy in polys.items():
        # gdstk method
        if len(xy) > 0:
            for p in gdstk.boolean(xy, clip_poly, operation, precision=1e-4):
                self._put_polygon(connect=(0,0,0),
                                  polygon=nazca.Polygon(points=p.points, layer=lay) )
        
        # # Grow and merge to heal small gaps
        # xy = [ nazca.clipper.\
        #     grow_polygons([poly], grow=1e-3, accuracy=1e-4, jointype='miter')[0]\
        #         for poly in xy ]
        # xy = nazca.clipper.merge_polygons(xy, accuracy=1e-4)
        
        # # Diff then add
        # for p in nazca.clipper.diff_polygons(xy, clip_poly, accuracy=1e-4):
        #     self._put_polygon(connect=(0,0,0),
        #                       polygon=nazca.Polygon(points=p, layer=lay) )
    
    return self
nazca.Cell.polydiff_layer = cell_polydiff_layer

# Invert layers of a cell according to existing layer or bounding box
def cell_invert(self, bounds, mod_layers, name=None, instantiate=None, keepother=False):
    """Build new cell with polygons in specified layer(s) inverted; implicitly flattens mod_layer.
    Also callable as `nazca.cell.invert`.
    
    Args:
        bounds (int, str, or polygon): layers or polygon to define inversion boundary
        mod_layer (int or str): layer(s) to invert; mod_layer' = bounds - mod_layer
        keepother (bool): keep unmodified layers
    
    Returns:
        Cell: new cell with clipped polygons
    """
    # Process options
    layers = layerlist(mod_layers)
    
    if name is None:
        name = self.cell_name+'_invert'+str(layers)
    
    if instantiate is None:
        instantiate = self.instantiate
    
    # Convert bounds to polygon if layers are specified
    if isinstance(bounds, int) or len(bounds[0]) == 1:
        bounds = layerlist(bounds)
        bounds = [self.get_polygons(l) for l in bounds]
    else:
        bounds = [bounds]   # Should be a list of polygons
    # Make sure 0-depth or 1-depth list of points
    while isinstance(bounds[0][0][0], list) or isinstance(bounds[0][0][0], tuple):
        bounds = bounds[0]
    
    # Merge bounds
    bounds = nazca.clipper.merge_polygons(bounds, accuracy=1e-8)
    
    # Get relevant polygons in layer(s)
    polys = nazca.defaultdict(list)
    for lay in layers:
        laypoly = self.get_polygons(lay)
        if len(laypoly) > 0:
            polys[lay] = laypoly
    if len(polys) < 1: # Nothing to do
        nazca.main_logger(
            msg='%s.invert(): No polygons found in mod_layers %s; returning %s.' % (self.cell_name, str(layers), self.cell_name), 
            level='warning')
        return self
    
    # Clip polygons
    for lay, xy in polys.items():
        # gdstk method
        # Clip
        polys[lay] = [p.points for p in 
                        gdstk.boolean(bounds, xy, 'not', precision=1e-4) ]
        # Merge results
        polys[lay] = [p.points for p in 
                        gdstk.boolean(polys[lay], polys[lay], 'or', precision=1e-4) ]
        
        # Make sure polygons are in reasonable order for Nazca's GDS writer
        # polys[lay] = nazca.clipper.merge_polygons(polys[lay], accuracy=1e-4)
    
    # New cell
    with nazca.Cell(name=name, instantiate=instantiate) as c:
        # Put unmodified layers if desired
        if keepother:
            self.remove_layer(layers).put('org', 0)
        
        # Add polygons
        for lay, xy in polys.items():
            # xy = flatten(xy)
            for poly in xy:
                nazca.Polygon(points=poly, layer=lay).put(0)
    
    # TODO: copy pins; this probably trips up pin connections?
    c.pin = self.pin
    c.default_in = self.default_in
    c.default_out = self.default_out
    
    # Copy other properties
    c.length_geo = self.length_geo if hasattr(self, 'length_geo') else None
    
    return c
nazca.Cell.invert = cell_invert

# Diff/clip layers of a cell
def cell_diff_layers(self, clip_layer=1, mod_layer=1001):
    """Modify current cell by clipping layers.
    Also callable as `nazca.cell.diff_layers`.
    
    Args:
        clip_layer (int or str): layer defining negative clip, will not change
        mod_layer (int or str): layer to modify; layer2' = layer2 - layer1
    
    Returns:
        Cell: modified self
    """
    # Get polygons to use as inverse clip
    clipPoly = self.get_polygons(clip_layer)
    if len(clipPoly) == 0:   # No modification
        return self
    
    # Perform difference
    self.polydiff_layer(clip_poly=clipPoly, mod_layer=mod_layer)
    
    return self
nazca.Cell.diff_layers = cell_diff_layers



### Custom bends
## Raised Sine s-curve
#   K. L. Kruse and C. T. Middlebrook, https://doi.org/10/ggp757
#   Calculation already implemented in nazca.generic_bend.sinebend_point
# Explicitly specified
def Tp_sinbend_dw(self, length=100, **kwargs):
    """Create a variable-width sine bend interconnect.
    Also callable as `nazca.Interconnect.sinebend_dw`.
    
    Args:
        length (float): interconnect length (dX)
        offset (float): interconnect offset (dY)
        width1 (float): begin width
        width2 (float): end width
        xs (str): xsection name
        layer (str | tuple | int): mask layer
        name (str): element name (default='viper')
        N (int): number of polygon points (default 1/µm)
        arrow (bool): draw connection arrows
        **kwarg: free parameters
    
    Returns:
        Cell: sine bend element
    """
    if length <= 0:
        return self.strt(length=0)
    
    def x(t, length, **kwargs): return t * length
    def y(t, offset, **kwargs): return offset * t - offset / (2 * pi) * sin(2 * pi * t)
    def w(t, width1, width2, **kwargs): 
        return width1 + (width2-width1) * t - (width2-width1) / (2 * pi) * sin(2 * pi * t)
    
    if 'N' not in kwargs:
        kwargs['N'] = max([round(length*5), 10])
    c = self.Tp_viper(x=x, y=y, w=w, 
                      params={'length': length, 'offset': 0})(**kwargs)
    
    # For some reason must explicitly set the instantiate flag; rebuild or with nazca.Cell()... doesn't work
    c.instantiate = False
    return c
nazca.interconnects.Interconnect.sinebend_dw = Tp_sinbend_dw

# Minimum length of a sine bend for a given offset and radius
def sinebend_minlen(dy, radius):
    """Minimum length of a sine bend for a given offset and radius."""
    return (2*pi*abs(dy)*abs(radius))**0.5


## Simple bend with zero-curvature endpoints
def ic_euler_bend(self, angle=180, radius=None, width=None, 
                  arrow=True, instantiate=True):
    """Simple Euler bend of total <angle> with zero-curvature endpoints.
    Also callable as `nazca.Interconnect.euler_bend`.
    
    Args:
        angle (float): Total angle to traverse
        radius (float): Minimum radius of curvature; default self.radius
    
    Returns:
        Cell: properly pinned waveguide cell with total waveguide length <length>
    """
    if radius is None:
        radius = self._getradius(nazca.cp.here(), radius, self.xs)
    
    # Name and deduplicate
    name = str(self.xs)+'.Euler.'+str([angle, radius])
    if name in nazca.cfg.cellnames.keys():
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(instantiate=False) as segment:
        if angle != 0:
            b = self.euler(radius=radius, angle=angle/2.0, width=width, arrow=False)
        else:
            b = self.strt(length=0, arrow=False)
        b.put(0).raise_pins(['a0'])
        b.put('b0', flip=True).raise_pins(namesin=['a0'], namesout=['b0'])
        
        segment.length_geo = 2*b.length_geo
        
        # Draw pins
        if arrow:
            for _, p in segment.ic_pins():
                nazca.make_pincell().put(p)
    
    return segment.flatten(name=name, instantiate=instantiate)
nazca.interconnects.Interconnect.euler_bend = ic_euler_bend


## Hairpin bend
def ic_hairpin(self, radius=None, overangle=50, direction='left', 
               arrow=True, instantiate=True):
    """Create a hairpin bend with continuous curvature via `euler` & `bend`.
    Also callable as `nazca.Interconnect.hairpin`.
    
    Args:
        radius (float): Minimum radius of curvature
        overangle (str): Extra angle above 180 degrees to hairpin
        direction (str): 'left' or 'right'
    
    Returns:
        Cell: properly pinned waveguide cell
        length = pi*radius*(1 + 4*overangle/180)
    """
    # Definitions
    edge = 180
    
    # Process inputs
    if radius is None:
        radius = self._getradius(nazca.cp.here(), radius, self.xs)
    if 'right' in direction:
        overangle *= -1
        edge *= -1
    
    # Name and deduplicate
    name = str(self.xs)+'.hairpin.'+str(direction)+'.'+str([abs(overangle), radius])
    if name in nazca.cfg.cellnames.keys():
        return nazca.cfg.cellnames[name]
    
    with nazca.Cell(instantiate=False) as segment:
        # Define bends
        b1 = self.euler_bend(angle=-overangle/2.0, radius=radius, arrow=False)
        b2 = self.euler(angle=overangle/2.0, radius=radius, arrow=False)
        b3 = self.bend(angle=edge, radius=radius, arrow=False)
        segment.length_geo = 2*b1.length_geo + 2*b2.length_geo + b3.length_geo
        
        # Place bends
        b1.put(0).raise_pins(['a0'])
        b2.put()
        b3.put()
        b = b2.put('b0', flip=True)
        b1.put(b.pinin).raise_pins(['b0'])
        
        # Draw pins
        if arrow:
            for _, p in segment.ic_pins():
                nazca.make_pincell().put(p)
    
    return segment.flatten(name=name, instantiate=instantiate)
nazca.interconnects.Interconnect.hairpin = ic_hairpin

