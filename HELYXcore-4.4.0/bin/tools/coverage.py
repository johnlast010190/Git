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
#    (c) 2024 Engys Ltd.

# Script
#
#
# Description
#     This script checks the coverage of features in HELYX
#
# ------------------------------------------------------------------------------


import os
import sys
import argparse

path_features = "../../etc/coverageLists/"

# Dictionary of features and description as defined in the coverage lists
TAG_STRINGS = {
    'bcGUI': 'Boundary conditions in HELYX-GUI',
    'bcTESTS': 'Boundary conditions in helyxTest',
    'buoyancyGUI': 'Buoyancy models in HELYX-GUI',
    'buoyancyTESTS': 'Buoyancy models in helyxTest',
    'motionGUI': 'Motion in HELYX-GUI',
    'motionTESTS': 'Motion in helyxTest',
    'turbGUI': 'Turbulence models in HELYX-GUI',
    'turbTESTS': 'Turbulence models in helyxTest',
    'statesGUI': 'States in HELYX-GUI',
    'statesTESTS': 'States in helyxTest',
    'statesModulesGUI': 'States with modules in HELYX-GUI',
    'statesModulesTESTS': 'States with modules in helyxTest',
}

# Reading features from files (temporary)
# These files will be probably created based on 
# what is in HELYX-CORE
def process_files(folder, tags):
    files_and_content = []
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        if os.path.isfile(file_path):
            try:
                with open(file_path, 'r') as file:
                    comment = file.readline().strip().lstrip('#').strip()
                    content = file.read()
                    words = content.split()
                    files_and_content.append((filename, comment, words))
            except FileNotFoundError:
                print(f"File not found: {file_path}")
            except Exception as e:
                print(f"An error occurred while processing {filename}: {e}")

    # Check if all specified tags are present in TAG_STRINGS
    invalid_tags = [tag for tag in tags if tag not in TAG_STRINGS]
    if invalid_tags:
        print(f"Error: Invalid tags specified: {', '.join(invalid_tags)}")
        print("Possible tags features:")
        for valid_tag in TAG_STRINGS:
            print(f"  {valid_tag}")
        sys.exit(1)

    tag_strings = [TAG_STRINGS.get(tag, tag) for tag in tags]

    # Filter files_and_content based on tags
    filtered_files_and_content = [
        file_info for file_info in files_and_content
        if any(tag.lower() in file_info[1].lower() for tag in tag_strings)
    ]

    return filtered_files_and_content


def print_possible_tags():
    print("Possible tags:")
    for valid_tag, description in TAG_STRINGS.items():
        print(f"  {valid_tag}")
    sys.exit(0)


def checkModels(folder, filename, strings_to_check):
    matching_files = []
    for root, dirs, files in os.walk(folder):
        if filename in files:
            file_path = os.path.join(root, filename)
            try:
                with open(file_path, 'r') as file:
                    content = file.read()
                    if any(word.lower() in content.lower() for word in strings_to_check):
                        matching_files.append(file_path)
            except Exception as e:
                print(f"An error occurred while processing {file_path}: {e}")
    return matching_files


def process_coverage(args):
    if args.tags is not None and not args.tags:
        print_possible_tags()
        sys.exit(0)
    else:
        if args.paths:
            info = process_files(path_features, args.tags or [])

            # Dictionary to store present/missing features
            missing_features = {}
            present_features = {}

            path_cases = {}

            for i in range(len(args.paths)):
                path = args.paths[i]
                print("\n")
                print(f"Features not present in: {path}:\n")
                if os.path.exists(path):
                    missing = []
                    present = []
                    paths = []
                    for feature in info[0][2]:
                        matching_files = checkModels(path, args.dict, [feature])
                        feature_printed = False  # Flag to track if the feature has been printed
                        for file_path in matching_files:
                            if feature not in present:  # Check if feature has already been printed
                                present.append(feature)
                                paths.append(file_path)
                                # print(f"{feature}")
                                feature_printed = True
                                break  # Exit the loop if feature is found once

                        if not feature_printed:
                            missing.append(feature)
                            print(f"{feature}")

                    missing_features[path] = missing
                    present_features[path] = present
                    path_cases[path] = paths
                else:
                    print(f"Invalid path specified: {path}")

            if len(args.paths) > 1:
                common_features = set(present_features[args.paths[0]])
                for path in args.paths[1:]:
                    common_features &= set(present_features.get(path, []))

                missing = [(index, feature) for index, feature in enumerate(present_features[args.paths[0]]) if feature not in common_features]
                print("\nFeatures present in " + args.paths[0] +
                    " and not in " + args.paths[1] + "\n")
                for index, feature in missing:
                    print(f"Feature: {feature}, Path: {os.path.dirname(os.path.dirname(path_cases[args.paths[0]][index]))}")
            else:
                print("\nOnly one repo specified, no comparison of coverage made.")
        else:
            parser.print_help()



def main():
    parser = argparse.ArgumentParser(
        description='Process files and search for tags.'
    )
    parser.add_argument(
        "-t",
        "--tags",
        nargs="*",
        help="Tags to features type"
    )
    parser.add_argument(
        "-p",
        "--paths",
        nargs="*",
        help="Paths to process"
    )
    parser.add_argument(
        "-d",
        "--dict",
        default="caseSetupDict",
        help="Specify the dictionary keyword to search for features (default: caseSetupDict)"
    )

    args = parser.parse_args()

    # Check coverage
    process_coverage(args)



if __name__ == '__main__':
    main()