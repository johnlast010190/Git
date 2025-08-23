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
    (c) 2021 OpenFOAM Foundation

Class
    Foam::systemDict

Description
    Read and return the specified dictionary from system or from path provided
    with the -dict option.

SourceFiles
    systemDict.C

\*---------------------------------------------------------------------------*/

#ifndef systemDict_H
#define systemDict_H

#include "global/argList/argList.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

IOobject systemDictIO
(
    const word& dictName,
    const argList& args,
    const objectRegistry& ob,
    const word& regionName = polyMesh::defaultRegion
);

IOdictionary systemDict
(
    const word& dictName,
    const argList& args,
    const objectRegistry& ob,
    const word& regionName = polyMesh::defaultRegion
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
