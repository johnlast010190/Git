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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::ensightParts

Description
    A collection of several ensightPart elements

SourceFiles
    ensightParts.C

\*---------------------------------------------------------------------------*/

#ifndef ensightParts_H
#define ensightParts_H

#include "containers/LinkedLists/user/SLPtrList.H"
#include "ensight/part/ensightPart.H"
#include "ensight/part/ensightPartFaces.H"
#include "ensight/part/ensightPartCells.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ensightParts Declaration
\*---------------------------------------------------------------------------*/

class ensightParts
:
    public SLPtrList<ensightPart>
{
public:

    //- Storage type used
    typedef SLPtrList<ensightPart> StorageType;


    // Constructors

        //- Construct from polyMesh
        ensightParts(const polyMesh&);

        //- Disallow default bitwise copy construction
        ensightParts(const ensightParts&) = delete;


    //- Destructor
    ~ensightParts();


    // Member functions

        //- Clear old information and construct anew from polyMesh
        void recalculate(const polyMesh&);

        //- Number of parts
        using StorageType::size;

        //- Write the geometry
        void write(ensightGeoFile&) const;

        //- Write summary information about the objects
        void writeSummary(Ostream&) const;

        //- Print various types of debugging information
        void dumpInfo(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ensightParts&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
