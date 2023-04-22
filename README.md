# photonics-pdk
Utilities and a PDK for a custom GaAs processe provided as Python packages for use in [nazca-degin](https://nazca-design.org/) ([fork with fixes](https://github.com/nickersonm/nazca-design)).


## Description

A utility library and PDK for a custom GaAs process, for use with the [nazca-degin](https://nazca-design.org/) GDS layout library. My [nazca fork with fixes](https://github.com/nickersonm/nazca-design) will be most compatible. KLayout scripts for postprocessing the nazca output are also included.


## Install

Clone this repository and add it to your Python PATH.

### Dependencies

- The [nazca-degin](https://nazca-design.org/) GDS layout library is the primary dependency. My [nazca fork with fixes](https://github.com/nickersonm/nazca-design) will be most compatible.
- numpy
- [gdstk](https://heitzmann.github.io/gdstk/reference_python.html) is used in `pdk_Fab6\pdk_05_functions.py` for additional speed, but can be replaced with nazca's `clipper` functions in the appropriate sections.
- [KLayout](https://www.klayout.de) is needed for some postprocessing.


## Usage

Sample mask definition files are provided in [`samples/`](./samples/). Make sure to set the KLayout path on the last line.

**TODO**

