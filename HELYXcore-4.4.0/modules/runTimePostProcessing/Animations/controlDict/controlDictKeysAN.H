/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : Dev
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
    (c) 2024 Engys Ltd.

Namespace
    Foam::functionObjects::animations::controlDictKeysAN

Description
    Contains the definition of the keys used in the controlDict for
    the animation FO

SourceFiles
    <controlDictKeysAN.C>

\*---------------------------------------------------------------------------*/

#pragma once

namespace Foam::functionObjects::animations::controlKeys
{
    extern const char* FRAME_RATE_KEY;
    extern const char* WIDTH_KEY;
    extern const char* HEIGHT_KEY;
    extern const char* OUTPUT_FORMAT_KEY;
    extern const char* TIMING_MODE_KEY;

    extern const char* TIME_START_KEY;
    extern const char* TIME_END_KEY;
    extern const char* TIME_SCALE_KEY;
    extern const char* FRAME_DURATION_KEY;
    extern const char* RTPP_INPUTS_KEY;
    extern const char* FFMPEG_PATH_KEY;
} // End namespace Foam
