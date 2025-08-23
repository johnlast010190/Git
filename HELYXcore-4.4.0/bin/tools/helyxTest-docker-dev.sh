#!/usr/bin/env bash

#---------------------------------------------------------------------------
#|       o        |
#|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
#|   o   O   o    |  Version : Dev
#|    o     o     |  ENGYS Ltd. <http://engys.com/>
#|       o        |
#---------------------------------------------------------------------------
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
#    (c) 2019-2024 Engys Ltd.
#
# Script
#     helyxTest-docker-dev.sh
#
# Description
#     Creates a standard docker environment for building HELYXcore.
#
#     The base image is an OpenSUSE v15.1. It should install everything needed
#     to build HELYXcore, along with helyxTest, which will be already in the
#     PATH. Note that cases are not downloaded and you'd need to use
#     `helyxTest init` to do so. If you have the cases already downloaded,
#     `init` again won't be needed, 'cause the user HOME directory is mounted
#     inside the Docker image, giving access to all the normal things.
#------------------------------------------------------------------------------

sshkey=""

print_help()
{
help_string=$(cat <<END
--------------------------- helyxTest-docker-dev help --------------------------

Creates a standard docker environment for building HELYXcore.

The base image is an OpenSUSE v15.1. It should install everything needed to
build HELYXcore, along with helyxTest and helyxVerify, both already in PATH.

Note that helyxTests are not downloaded and you'd need to use 'helyxTest init'
to do so. If you have the cases already downloaded, 'init' again won't be
needed, 'cause the user HOME directory is mounted inside the Docker image,
giving access to all the normal things.

If you have any issues using helyxTest to download the cases, try doing it so
outside the image.

-----------------------------------  Usage  ------------------------------------

    helyxTest-docker-dev [<flags>]

    flags
        All flags are optional, and the available flags are given below.
        All flags support both camelCase (like OpenFOAM) and hyphenation.

  -a | --auth
      Authenticate the machine on the Docker registry. Needs to be
      used only once.

        -h | --help
            Shows this text.

--------------------------------------------------------------------------------
END
)
echo "$help_string"
exit 0
}

OPTS=$(getopt -a \
-o "h:,a,c" \
-l "help:,auth:,clean" \
-n 'helyxTest-docker-dev' -- "$@")

# shellcheck disable=SC2181
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

eval set -- "$OPTS"

while true; do
    case "$1" in
      -h | --help )
          print_help
          ;;
      -a | --auth )
        docker login registry.engys.dev/ops -u ops -p $TOKEN
        shift
        ;;
      -c | --clean )
        echo "Removing old, unused developer images..."
        old_images=$(docker image ls registry.engys.dev/ops/jenkins-images/evtk --filter dangling=true -q)
        if [ ! -z "$old_images" ]; then
          docker image rmi $old_images > /dev/null
        fi
        echo "Images removed."
        exit 0
        ;;
      -- )
          if [ $# -gt 1 ]; then
              print_help
          fi
          shift
          break
          ;;
      * )
          break ;;
    esac
done

if ! command -v docker &> /dev/null
then
  echo "Docker is required."
  exit 1
fi

image_name="registry.engys.dev/ops/jenkins-images/evtk:latest"
docker pull "${image_name}" &&
  docker run -it \
    -u $(id -u):$(id -g) \
    -e HOME=${HOME} \
    -v /etc/group:/etc/group:ro \
    -v /etc/passwd:/etc/passwd:ro \
    -v /etc/shadow:/etc/shadow:ro \
    -v /etc/sudoers.d:/etc/sudoers.d:ro \
    -v ${HOME}:${HOME}:rw \
    -w ${HOME} \
    --security-opt seccomp=unconfined \
    --rm \
    ${image_name} \
    /bin/bash
