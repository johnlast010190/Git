/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.0.1
|    o     o     |  ENGYS Ltd. <http://engys.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of HELYXcore.
    HELYXcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    HELYXcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HELYXcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HELYXcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2022-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "colours.H"

namespace Foam::functionObjects::runTimeVis
{

const scalar ENGYS_CYAN_DARK[3] = { static_cast<scalar>(115 / 255.0), static_cast<scalar>(160 / 255.0), static_cast<scalar>(200 / 255.0) };
const scalar ENGYS_BLUE_LIGHT[3] = { static_cast<scalar>(40 / 255.0), static_cast<scalar>(78 / 255.0), static_cast<scalar>(117 / 255.0) };
const scalar BLACK[3] = { 0, 0, 0 };
const scalar GRAY[3] = { static_cast<scalar>(150 / 255.0), static_cast<scalar>(150 / 255.0), static_cast<scalar>(150 / 255.0) };
const scalar MAGENTA[3] = { 1, 0, 1 };

const scalar NAN_BAR_COLOR[3] = {ENGYS_BLUE_LIGHT[0], ENGYS_BLUE_LIGHT[1], ENGYS_BLUE_LIGHT[2]};
const scalar NAN_ACTOR_COLOR[3] = {ENGYS_CYAN_DARK[0], ENGYS_CYAN_DARK[1], ENGYS_CYAN_DARK[2]};

} // End namespace
