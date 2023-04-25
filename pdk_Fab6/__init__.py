#!/usr/bin/env python3
# 
# Provides components for GaAs process "Fab6"
# 
# @Authors: Michael Nickerson
# @email: nickersonm@ece.ucsb.edu
# 2023(c)
# 

"""
GaAs components with 'Fab6' process
(c) Michael Nickerson 2023
"""

# Allow import of sibling packages
from sys import path
from os.path import abspath
path.insert(0, abspath('..'))

# Make sure nazca is loaded
import nazca

# Contained modules
from .pdk_20_technology import *
from .pdk_30_components import *
from .pdk_35_devices import *
from .pdk_40_devutil import *
