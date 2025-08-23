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
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::polyBoundaryMeshEntries

Description
    Foam::polyBoundaryMeshEntries

SourceFiles
    polyBoundaryMeshEntries.C

\*---------------------------------------------------------------------------*/

#ifndef polyBoundaryMeshEntries_H
#define polyBoundaryMeshEntries_H

#include "db/regIOobject/regIOobject.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "db/dictionary/entry/entry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyBoundaryMeshEntries
:
    public regIOobject,
    public PtrList<entry>
{

public:

    //- Runtime type information
    TypeName("polyBoundaryMesh");


    // Constructors

        explicit polyBoundaryMeshEntries(const IOobject& io)
        :
            regIOobject(io),
            PtrList<entry>(readStream(typeName))
        {
            close();
        }


   // Member functions

        bool writeData(Ostream&) const
        {
            NotImplemented;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
