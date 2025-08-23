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
    (c) 2011 OpenFOAM Foundation

InClass
    Foam::writeFunctions

Description

SourceFiles
    writeFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef writeFunctions_H
#define writeFunctions_H

#include "primitives/Scalar/scalar/scalar.H"
#include "primitives/strings/string/string.H"
#include "primitives/ints/lists/labelList.H"
#include <iostream>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool writeStr80(std::ostream& os, const string& str);
bool writeInt(std::ostream& os, int val1);
bool writeFloat(std::ostream& os, scalar val1);
void writeBytes(char* start, int nBytes);
void writeWallFlags(Ostream& os, label cellI, const labelList& wallFlags);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
