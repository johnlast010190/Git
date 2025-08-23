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
#
# Description
#     Helper classes for getting information about the HELYX source tree.
#
# ------------------------------------------------------------------------------

import pathlib


class FilePathGetter:
    def __init__(self, parent_directory, subdirs_to_search):
        self.parent_dir = parent_directory
        self.subdirs_to_search = []
        for subdir_name in subdirs_to_search:
            subdir = pathlib.Path(parent_directory / subdir_name)
            if subdir.exists():
                self.subdirs_to_search.append(subdir)
            else:
                print(
                    f"Warning:  Could not find subdirectory {subdir_name} on path {parent_directory}"
                )

        self.all_file_paths = []
        # First, look for top-level CMakeLists.txt file
        if (self.parent_dir / "CMakeLists.txt").is_file():
            self.all_file_paths.append(self.parent_dir / "CMakeLists.txt")
        # Now iterate through for all other files
        for subdir in self.subdirs_to_search:
            temp_list_of_files = [x for x in list(subdir.rglob("*")) if x.is_file()]
            self.all_file_paths += temp_list_of_files


def replace_string_in_file(input_file, string_to_replace, replacement):
    with open(input_file, "r") as f:
        content = f.read()

    content = content.replace(string_to_replace, replacement)
    with open(input_file, "w") as f:
        f.write(content)


class HelyxDirectoryInfo:
    def __init__(self):
        # Assume this file is kept in bin/tools, derive helyx_project_dir.  Use
        # the same names as the environment variables (even though they're a bit
        # inconsistent)
        self.helyx_project_dir = pathlib.Path(__file__).absolute().parents[2]
        self.helyx_modules = self.helyx_project_dir / "modules"
        self.adjoint_dir = self.helyx_modules / "HELYX-adjoint"
        self.coupled_dir = self.helyx_modules / "HELYX-coupled"
        self.hydro_dir = self.helyx_modules / "HELYX-hydro"
        self.marine_dir = self.helyx_modules / "HELYX-marine"
        self.waves2Foam_dir = self.helyx_modules / "waves2Foam"
        self.runTimePostProcessing_dir = self.helyx_modules / "runTimePostProcessing"
        #self.examples = self.helyx_project_dir / "examples"
        #self.unitTests = self.helyx_project_dir / "unitTests"

        # Suspect messing with swak is probaly a bad idea
        # swak4Foam_dir = modules_dir / 'swak4Foam'

        # This is inelegant...
        self.directories_dict = {
            self.helyx_project_dir: [
                "applications",
                "bin",
                "doc",
                "etc",
                "src",
                "thirdParty",
                "tutorials",
                "examples",
                "unitTests"
            ],
            self.adjoint_dir: [
                "applications",
                "bin",
                "src",
                "tutorials",
            ],
            self.coupled_dir: [
                "applications",
                "src",
                "tutorials",
            ],
            self.hydro_dir: [
                "applications",
                "src",
                "tutorials",
            ],
            self.marine_dir: [
                "applications",
                "src",
                "tutorials",
            ],
            self.waves2Foam_dir: [
                "applications",
                "bin",
                "doc",
                "src",
                "tutorials",
            ],
            self.runTimePostProcessing_dir: [
                "FieldSampling",
                "PostDictObjectProviderDatabase",
                "RunTimeVisualisation",
                "SurfaceStatistics",
            ],
        }


if __name__ == "__main__":
    directory_info = HelyxDirectoryInfo()
    for parent_dir in directory_info.directories_dict:
        print(f"Top-level directory {parent_dir}:")
        print(f"    {directory_info.directories_dict[parent_dir]}:")
        # test = FilePathGetter(parent_dir, directory_info.directories_dict[parent_dir])
