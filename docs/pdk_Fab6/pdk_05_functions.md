Module pdk_Fab6.pdk_05_functions
================================
Utility functions for pdk_Fab6.
(c) Michael Nickerson 2023

Functions
---------

    
`Tp_sinbend_dw(self, length=100, **kwargs)`
:   Create a variable-width sine bend interconnect.
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

    
`cellArray(cells, space=20, nx=3)`
:   Generic array of cells.

    
`cellChangeName(cell, newname)`
:   Change a cell name.

    
`cell_diff_layers(self, clip_layer=1, mod_layer=1001)`
:   Modify current cell by clipping layers.
    Also callable as `nazca.cell.diff_layers`.
    
    Args:
        clip_layer (int or str): layer defining negative clip, will not change
        mod_layer (int or str): layer to modify; layer2' = layer2 - layer1
    
    Returns:
        Cell: modified self

    
`cell_flatten(self, name=None, instantiate=None)`
:   Rebuild cell to merge all polygons per layer.
    Also callable as `nazca.cell.flatten`.
    
    Args:
        name (str): new name for the cell
        instantiate (bool): instantiate the new cell?
    
    Returns:
        Cell: new cell with merged polygons

    
`cell_get_polygons(self, layers)`
:   Return all polygons in <layers> as point lists.
    Also callable as `nazca.cell.get_polygons`.
    
    Args:
        layers (list of str): names of layers to retreive polygons from
    
    Returns:
        List: list of polygon points

    
`cell_invert(self, bounds, mod_layers, name=None, instantiate=None, keepother=False)`
:   Build new cell with polygons in specified layer(s) inverted; implicitly flattens mod_layer.
    Also callable as `nazca.cell.invert`.
    
    Args:
        bounds (int, str, or polygon): layers or polygon to define inversion boundary
        mod_layer (int or str): layer(s) to invert; mod_layer' = bounds - mod_layer
        keepother (bool): keep unmodified layers
    
    Returns:
        Cell: new cell with clipped polygons

    
`cell_polydiff_layer(self, clip_poly, mod_layer=1001, operation='not')`
:   Rebuild cell to clip all polygons in specified layer(s); implicitly flattens mod_layer.
    Also callable as `nazca.cell.polydiff_layer`.
    
    Args:
        clip_poly (list): polygon geometry to use as clip
        mod_layer (int or str): layer(s) to modify; mod_layer' = mod_layer - clip_poly
    
    Returns:
        Cell: existing cell with clipped polygons

    
`cell_remove_layer(self, layers, name=None, instantiate=None)`
:   Remove all polygons in `layers`.
    Also callable as `nazca.cell.remove_layer`.
    
    Args:
        layers (list of str or int or xs): names of layers to delete polygons from
        name (str): new name for the cell
        instantiate (bool): instantiate the new cell?
    
    Returns:
        Cell: new cell with removed polygons

    
`ic_euler_bend(self, angle=180, radius=None, width=None, arrow=True, instantiate=True)`
:   Simple Euler bend of total <angle> with zero-curvature endpoints.
    Also callable as `nazca.Interconnect.euler_bend`.
    
    Args:
        angle (float): Total angle to traverse
        radius (float): Minimum radius of curvature; default self.radius
    
    Returns:
        Cell: properly pinned waveguide cell with total waveguide length <length>

    
`ic_hairpin(self, radius=None, overangle=50, direction='left', arrow=True, instantiate=True)`
:   Create a hairpin bend with continuous curvature via `euler` & `bend`.
    Also callable as `nazca.Interconnect.hairpin`.
    
    Args:
        radius (float): Minimum radius of curvature
        overangle (str): Extra angle above 180 degrees to hairpin
        direction (str): 'left' or 'right'
    
    Returns:
        Cell: properly pinned waveguide cell
        length = pi*radius*(1 + 4*overangle/180)

    
`parsePin(pin, default='out', rot=False)`
:   Parse a pin/location/node into a real Pin object, with all transformations applied.

    
`pinGeoDiff(p1, p2)`
:   Get the geometric length between two pins with `nazca.cell.geolen`.

    
`plotPolygon(p)`
:   Simple debugging function to plot polygons.

    
`sinebend_minlen(dy, radius)`
:   Minimum length of a sine bend for a given offset and radius.

    
`taperLen(w1, w2)`
: