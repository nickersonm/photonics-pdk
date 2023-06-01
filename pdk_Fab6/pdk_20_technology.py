#!/usr/bin/env python3
# 
# Provides waveguide cross-sections for Fab6 process
# 
# @Authors: Michael Nickerson
# @email: nickersonm@ece.ucsb.edu
# 2023(c)
# 

"""
Waveguide cross-sections with 'Fab6' process.
(c) Michael Nickerson 2023
"""

import nazca



### Define layers, i.e. fab process steps
nazca.add_layer(name='ProtectRidge', layer=1, accuracy=1e-3, \
    remark='Protect against deep ridge etch')

nazca.add_layer(name='ProtectRegrowth', layer=15, accuracy=1e-3, \
    remark='Protect active regions against MQW selective etch')

nazca.add_layer(name='EtchRib', layer=2, accuracy=1e-3, \
    remark='Shallow rib etch regions')

nazca.add_layer(name='EtchIso', layer=3, accuracy=1e-3, \
    remark='Etch isolation-ridge regions')

nazca.add_layer(name='MetalVia', layer=4, accuracy=1e-3, \
    remark='Open P vias regions')

nazca.add_layer(name='MetalTop', layer=5, accuracy=1e-2, \
    remark='P-metal, top metal, and pad area')

nazca.add_layer(name='CleaveLane', layer=10, accuracy=1e-2, \
    remark='Remove oxide to allow easier cleave')

nazca.add_layer(name='CleaveKeepaway', layer=12, accuracy=1e-2, \
    remark='Protect from cleave lane')

nazca.add_layer(name='MetalKeepaway', layer=13, accuracy=1e-2, \
    remark='Temporary metal routing exclusion layer, to be differenced with MetalTop')

nazca.add_layer(name='AlignmentMark', layer=14, accuracy=1e-2, \
    remark='Alignment marks only')

# Informational layers
nazca.add_layer(name='DieArea', layer=1001, accuracy=1e-2, \
    remark='Die Area')
nazca.add_layer(name='DieReserved', layer=1005, accuracy=1e-2, \
    remark='Die Reserved')
nazca.add_layer(name='Annotation', layer=1111, accuracy=1e-1, \
    remark='Annotation layer')


### Cross sections
# Common sizes
textHeight = 35     # Minimum text height
keepawaySize = 50   # Cleave lane keepaway
contactWidth = 30   # Typical metal contact width
traceWidth = 20     # Typical metal trace width
metalBuffer = 15    # Metal-to-via buffer
viaInsetRib = 0.5   # Inset for top vias
viaInsetRidge = 0.75# Inset for top vias
typicalBuffer = 4   # Typical buffer size
cleaveWidth = 20    # Cleave lane width

# # Effective index for this epitaxy at 1030 nm, active and passive regions
# neff_A = lambda w : 0.0007082 * w + 3.4392
# neff_P = lambda w : -0.0379/(w**2) + 3.4410


## Passive sections
# Passive Deep WG
wgPassive = nazca.add_xsection(name='wgPassive', description='Deep ridge WG')
wgPassive.width = 2.0
wgPassive.radius = 125
wgPassive.minimum_radius = 125
wgPassive.taper = 20
nazca.add_layer2xsection(xsection='wgPassive', layer='ProtectRidge')
nazca.add_layer2xsection(xsection='wgPassive', layer='CleaveKeepaway', growx=keepawaySize, growy=keepawaySize)

# Passive Deep WG for bends only
wgBend = nazca.add_xsection(name='wgBend', description='Deep ridge WG')
wgBend.width = 1.5
wgBend.radius = 75
wgBend.minimum_radius = 75
wgBend.taper = wgPassive.taper
nazca.add_layer2xsection(xsection='wgBend', layer='ProtectRidge')
nazca.add_layer2xsection(xsection='wgBend', layer='CleaveKeepaway', growx=keepawaySize, growy=keepawaySize)

# Isolated WG
wgIso = nazca.add_xsection(name='wgIso', description='Isolated deep ridge WG')
wgIso.width = wgPassive.width
wgIso.radius = wgPassive.radius
wgIso.minimum_radius = wgPassive.minimum_radius
wgIso.taper = wgPassive.taper
wgIso.outset = 3.0  # SSA-type exposure
nazca.add_layer2xsection(xsection='wgIso', layer='ProtectRidge')
nazca.add_layer2xsection(xsection='wgIso', layer='EtchIso', growx=wgIso.outset)
nazca.add_layer2xsection(xsection='wgIso', layer='CleaveKeepaway', growx=keepawaySize, growy=keepawaySize)

# Passive Shallow WG
wgShallow = nazca.add_xsection(name='wgShallow', description='Shallow rib WG')
wgShallow.width = 3.0
wgShallow.pedestal = contactWidth + typicalBuffer
wgShallow.radius = 1e3
wgShallow.minimum_radius = 1e3
wgShallow.taper = 50
nazca.add_layer2xsection(xsection='wgShallow', layer='ProtectRidge', growx=(wgShallow.pedestal - wgShallow.width)/2, growy=0)
nazca.add_layer2xsection(xsection='wgShallow', layer='EtchRib', 
                         leftedge=( 0, (wgShallow.pedestal + typicalBuffer)/2 ), rightedge=( 0.5, 0 ))
nazca.add_layer2xsection(xsection='wgShallow', layer='EtchRib', 
                         rightedge=( 0, -(wgShallow.pedestal + typicalBuffer)/2 ), leftedge=( -0.5, 0 ))
nazca.add_layer2xsection(xsection='wgShallow', layer='CleaveKeepaway', growx=keepawaySize, growy=keepawaySize)

# Cleave lane
xsCleaveLane = nazca.add_xsection(name='xsCleaveLane', description='Cleave trench')
xsCleaveLane.width = cleaveWidth
nazca.add_layer2xsection(xsection='xsCleaveLane', layer='CleaveLane')
nazca.add_layer2xsection(xsection='xsCleaveLane', layer='MetalKeepaway', growx=metalBuffer, growy=metalBuffer)


## Phase Modulator
# Modulator: passive + metalP
wgModulator = nazca.add_xsection(name='wgModulator', description='Modulator waveguide portion')
wgModulator.width = 3
wgModulator.radius = wgPassive.radius
wgModulator.minimum_radius = wgPassive.minimum_radius
wgModulator.taper = wgPassive.taper
nazca.add_layer2xsection(xsection='wgModulator', layer='ProtectRidge')
nazca.add_layer2xsection(xsection='wgModulator', layer='CleaveKeepaway', growx=keepawaySize, growy=keepawaySize)
nazca.add_layer2xsection(xsection='wgModulator', layer='MetalVia', growx=-viaInsetRidge)
nazca.add_layer2xsection(xsection='wgModulator', layer='MetalTop', 
                         leftedge=( 0, contactWidth/2 ), rightedge=( 0, -contactWidth/2 ))


## Active section
# Laser or SOA: shallow + active + metalP
wgActive = nazca.add_xsection(name='wgActive', description='Active shallow rib waveguide')
wgActive.width = wgShallow.width
wgActive.pedestal = wgShallow.pedestal
wgActive.radius = 1e3
wgActive.minimum_radius = 1e3
wgActive.taper = 50
nazca.add_layer2xsection(xsection='wgActive', layer='ProtectRidge', growx=(wgActive.pedestal - wgActive.width)/2, growy=0)
nazca.add_layer2xsection(xsection='wgActive', layer='ProtectRegrowth', growx=(wgActive.pedestal - wgActive.width + 2*typicalBuffer)/2, growy=0)
nazca.add_layer2xsection(xsection='wgActive', layer='EtchRib', 
                         leftedge=( 0, (wgShallow.pedestal + typicalBuffer)/2 ), rightedge=( 0.5, 0 ))
nazca.add_layer2xsection(xsection='wgActive', layer='EtchRib', 
                         rightedge=( 0, -(wgShallow.pedestal + typicalBuffer)/2 ), leftedge=( -0.5, 0 ))
nazca.add_layer2xsection(xsection='wgActive', layer='CleaveKeepaway', growx=keepawaySize + (wgActive.pedestal - wgActive.width)/2, growy=keepawaySize)
nazca.add_layer2xsection(xsection='wgActive', layer='MetalVia', growx=-viaInsetRib)
nazca.add_layer2xsection(xsection='wgActive', layer='MetalTop', 
                         leftedge=( 0, contactWidth/2 ), rightedge=( 0, -contactWidth/2 ))

# Deep etch active region: passive + active + metalP
wgActiveDeep = nazca.add_xsection(name='wgActiveDeep', description='Active deep ridge waveguide')
wgActiveDeep.width = wgShallow.width
wgActiveDeep.radius = 150
wgActiveDeep.minimum_radius = 150
wgActiveDeep.taper = wgPassive.taper
nazca.add_layer2xsection(xsection='wgActiveDeep', layer='ProtectRidge')
nazca.add_layer2xsection(xsection='wgActiveDeep', layer='ProtectRegrowth', growx=typicalBuffer, growy=0)
nazca.add_layer2xsection(xsection='wgActiveDeep', layer='CleaveKeepaway', growx=keepawaySize, growy=keepawaySize)
nazca.add_layer2xsection(xsection='wgActiveDeep', layer='MetalVia', growx=-viaInsetRidge)
nazca.add_layer2xsection(xsection='wgActiveDeep', layer='MetalTop', 
                         leftedge=( 0, contactWidth/2 ), rightedge=( 0, -contactWidth/2 ))


## Metal regions
# Via plus top metal
metalP = nazca.add_xsection(name='metalP', description='P metal for top-P connections; includes via')
metalP.width = wgModulator.width
metalP.radius = wgModulator.width
metalP.drc_angle = False
metalP.drc_width = False
nazca.add_layer2xsection(xsection='metalP', layer='MetalVia', growx=-viaInsetRib)
nazca.add_layer2xsection(xsection='metalP', layer='MetalTop', 
                         leftedge=( 0, contactWidth/2 ), rightedge=( 0, -contactWidth/2 ))

# Generic top metal pad area
metalTop = nazca.add_xsection(name='metalTop', description='Top metal for routing, bonding, and probing')
metalTop.width = traceWidth
metalTop.radius = traceWidth
metalTop.drc_angle = False
metalTop.drc_width = False
nazca.add_layer2xsection(xsection='metalTop', layer='MetalTop')



### DRC for connections
nazca.cfg.drc_rule_xs = {
    'wgPassive': ['wgPassive', 'wgIso', 'wgModulator', 'wgBend'],
    'wgBend': ['wgPassive', 'wgBend'],
    'wgShallow': ['wgActive', 'wgShallow'],
    'wgModulator': ['wgModulator', 'wgIso', 'metalTop', 'metalP'],
    'wgActive': ['wgActive', 'wgShallow', 'metalTop', 'metalTop', 'metalP'],
    'wgActiveDeep': ['wgActiveDeep', 'metalTop', 'metalTop', 'metalP'],
    'metalTop': ['metalTop', 'metalP', 'wgModulator', 'wgActive', 'wgActiveDeep'],
    'metalP': ['metalTop', 'metalP', 'wgModulator', 'wgActive', 'wgActiveDeep']
}

