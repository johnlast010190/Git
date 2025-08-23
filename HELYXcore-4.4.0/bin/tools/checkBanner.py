#!/usr/bin/python3
# ---------------------------------------------------------------------------
# |       o        |
# |    o     o     |  HELYX (R) : Open-source CFD for Enterprise
# |   o   O   o    |  Version : 4.4.0
# |    o     o     |  ENGYS Ltd. <http://engys.com/>
# |       o        |
# ---------------------------------------------------------------------------
# License
#    This file is part of HELYXcore.
#    HELYXcore is based on OpenFOAM (R) <http://www.openfoam.org/>.
#
#    HELYXcore is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    HELYXcore is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#    for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with HELYXcore.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright
#    (c) 2021 Engys Ltd.

# Script
#
#
# Description
#     This script replaces an OpenFOAM banner with an HELYX one.
#
# ------------------------------------------------------------------------------

import pathlib
import re
from os import walk
from os.path import isfile, join

# check recursively files in this directory
path = str(pathlib.Path.home()) + "/HELYXcore-dev/tutorials/"
#path = str(pathlib.Path.home()) + "/HELYXcore-dev/etc/caseDicts"

version = "3.3.2"

HELYX_BANNER = [
f"/*--------------------------------*- C++ -*----------------------------------*\\\n",
f"|       o        |                                                            |\n",
f"|    o     o     |  HELYX (R) : Open-source CFD for Enterprise                |\n",
f"|   o   O   o    |  Version : {version}                                           |\n",
f"|    o     o     |  ENGYS Ltd. <http://engys.com/>                            |\n",
f"|       o        |                                                            |\n",
f"\\*---------------------------------------------------------------------------*/\n"
]

def check_banner(filename):
    # read only first 7 lines
    try:
        with open(filename) as f:
            banner = [next(f) for x in range(7)]
    #except (StopIteration):
    except:
        print(f"cannot read {filename}")
        return

    # if "F ield" and "O peration" found it's an old OpenFOAM banner
    foam_banner = False
    if (re.search("F ield",     banner[2]) and
        re.search("O peration", banner[3]) and
        re.search("\\*-------", banner[6])):
        foam_banner = True

    if foam_banner:
        print(f"Found OpenFOAM banner! I'll replace it with the ENGYS one")
        with open(filename) as f:
            lines = f.readlines()
        for n in range(7):
            lines[n] = HELYX_BANNER[n]
        with open(filename, "w") as f:
            f.writelines(lines)

# find files recursively
onlyfiles = list()
for (dirpath, dirnames, filenames) in walk(path):
    onlyfiles += [join(dirpath, file) for file in filenames]

for f in onlyfiles:
    filename = join(path, f)
    print(f"\nChecking for banner in {filename}")
    check_banner(filename)

