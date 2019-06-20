# ipkiss_meep_integration
![GitHub](https://img.shields.io/github/license/luceda/ipkiss_meep_integration.svg)

Integration between IPKISS (Luceda Photonics) and the Meep FDTD solver

IPKISS: www.lucedaphotonics.com
Meep: https://meep.readthedocs.io/en/latest/

Licensing
---------

* IPKISS is a commercially licensed tool by Luceda Photonics.
* Meep is licensed under GPLv2 
* This separate package is licensed under the MIT License, see LICENSE. 
* A valid license of IPKISS and an installation of Meep will be required in order to run it. Contact sales@lucedaphotonics.com for support.

Usage
-----
This package exports meep python scripts which can be run in Meep.

Warning: do not import meep directly within the code in this package, since that is incompatible with Meep's GPL license.
