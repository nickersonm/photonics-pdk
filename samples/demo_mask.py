#!/usr/bin/env python3
# 
# Demo mask for pdk_Fab6
# 
# @Authors: Michael Nickerson
# @url: http://nickersondevices.com
# 2023(c)
# 

import nazca
import pdk_Fab6 as PDK

with nazca.Cell(name='DFB') as cellDFB:
    PDK.inputLabel(text='PDK DFB, 2 mm', xs=None).put(50,0) # xs=None to avoid putting an input waveguide
    PDK.deviceDFB(length=2e3).put()
    PDK.soaToOut(x=5e3, soalen=500, angle=7).put()   # Straight waveguide to output at 5 mm border, with 500 µm SOA before 7° angled output

dieFile = 'demo_mask.gds'
nazca.export_gds(topcells=cellDFB, filename=dieFile, clear=False)

# Postprocess with KLayout
print('Postprocessing with KLayout...')
from subprocess import call
from os.path import abspath
postDRC = [abspath('..\\pdk_Fab6\\Fab6_postprocess.lydrc'),
           abspath('..\\pdk_Fab6\\Fab6_postprocess_merge.lydrc')]
for drc in postDRC:
    call(['klayout_app', '-b', '-r', drc, '-rd', 'input='+abspath(dieFile)])
