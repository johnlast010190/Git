#!/bin/bash

set -o errexit

# Historically, the FOAM buildsystem used "lnIncludes" during compilation.
# The build system has since moved on, but this script may be used to generate them in
# the source tree anyway, for those who want this.

if [[ $# -gt 1 ]]; then
  echo "-- Skiping lnInclude generation"
  echo "Provide only one path to generate the lnInclude folder, or leave it blank for HELYXcore-<version>"
  exit 0
elif [[ $# -eq 1 ]] && [ ! -d "$1" ]; then
  echo "-- Skiping lnInclude generation"
  echo "The path '$1' could not be found"
  exit 0
elif [[ $# -eq 1 ]] && [ -L "$1" ]; then
  echo "-- Skiping lnInclude generation"
  echo "The path '$1' should not be a symlink"
  exit 0
fi


function mkLnDir() {
  BASEDIR="${1}"
  LNDIR="${BASEDIR}/lnInclude"
  rm -Rf "$LNDIR"
  mkdir "$LNDIR"
  find "$BASEDIR" -type f \( -iname \*.cpp -o -iname \*.h -o -iname \*.hpp -o -iname \*.c -o -iname \*.type \) -print0 | xargs -0 realpath -z --relative-to "$LNDIR" | xargs -0 -I{} ln -sf {} "$LNDIR"
}


# For given path outside HELYXcore-<version>
if [[ $# -eq 1 ]]; then

  # mkLnDir "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
  #echo "-- Generating lnIncludes on: '$(realpath "${1}")'"
  mkLnDir "$(realpath "${1}")"

# For HELYXcore-<version>/
else

  #echo "-- Generating lnIncludes on: '$(dirname "$0")/../'"
  cd "$(dirname "$0")"/../

  # Identify where lnInclude directories should go. Some are hardcoded, the rest are governed by every
  # location that calls add_helyx_library()
  for i in $(find . -name CMakeLists.txt); do
    if grep -q add_helyx_library $i; then
      mkLnDir "$(dirname $i)"
    fi
  done

  # Some random directories were hardcoded. These correspond to those which are in need of refactoring anyway
  # (since they ought to be using `add_helyx_library`.
  mkLnDir src/OSspecific
  mkLnDir src/OpenFOAM
  mkLnDir src/TurbulenceModels/phaseCompressible
  mkLnDir src/TurbulenceModels/phaseIncompressible
  mkLnDir src/parallel/decompose/ptscotchDecomp

  cd -

fi
