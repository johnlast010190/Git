/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
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
    (c) 2007 Takuya OSHIMA <oshima@eng.niigata-u.ac.jp>

Description
    Compatibility symbols across different versions of OpenFOAM.

\*---------------------------------------------------------------------------*/

#ifndef convectiveOutletConfig_H
#define convectiveOutletConfig_H

// OF >= 1.4 requires a copy ctor
#if HELYX_PROJECT_VERSION_MAJOR > 1 || (HELYX_PROJECT_VERSION_MAJOR == 1 \
    && HELYX_PROJECT_VERSION_MINOR >= 4)
#define NEEDS_COPY_CTOR 1
#else
#define NEEDS_COPY_CTOR 0
#endif

// OF >= 1.4.1 requires ctor with DimensionedField arguments
#if HELYX_PROJECT_VERSION_MAJOR == 1 && (HELYX_PROJECT_VERSION_MINOR < 4 \
    || HELYX_PROJECT_VERSION_MINOR == 4 && HELYX_PROJECT_VERSION_PATCH == 0)
#define FIELD_IS_DIMENSIONED 0
#else
#define FIELD_IS_DIMENSIONED 1
#endif

#endif
