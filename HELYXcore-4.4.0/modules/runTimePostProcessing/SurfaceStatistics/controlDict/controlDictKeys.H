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
    (c) 2019 OpenCFD Ltd.
    (c) 2023 Engys Ltd.

Namespace
    Foam::functionObjects::surfaceStat::controlDictKeys

Description
    Contains the definition of the keys used in the controlDict for
    the surface statistics FO

SourceFiles
    <controlDictKeys.C>

\*---------------------------------------------------------------------------*/

#pragma once

namespace Foam::functionObjects::surfaceStat::controlKeys
{
    extern const char* FLIP_NORMAL_KEY;
    extern const char* WEIGHTING_KEY;
    extern const char* FIELDS_KEY;
    extern const char* SOURCES_KEY;
    extern const char* ORIENTATION_POINT_KEY;

    extern const char* EXPORT_RESULTS_KEY;
    extern const char* REFERENCE_FRAME_KEY;
} // End namespace Foam
