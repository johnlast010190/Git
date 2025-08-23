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
    (c) 2016 OpenCFD Ltd.
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::ensightPartFaces

Description
    An implementation of ensightPart to hold mesh faces.

SourceFiles
    ensightPartFaces.C

\*---------------------------------------------------------------------------*/

#ifndef ensightPartFaces_H
#define ensightPartFaces_H

#include "ensight/part/ensightPart.H"
#include "ensight/part/ensightFaces.H"
#include "db/typeInfo/typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ensightPartFaces Declaration
\*---------------------------------------------------------------------------*/

class ensightPartFaces
:
    public ensightFaces,
    public ensightPart
{
    // Private data

        //- Start offset for patch
        const label start_;

        //- Patch index
        const label patchIndex_;

        //- Faces referenced
        const faceList& faces_;

        //- pointField referenced
        const pointField& points_;

        //- Can skip local point renumbering when points are contiguous
        const bool contiguousPoints_;


    // Private Member Functions

        //- Track points used
        localPoints calcLocalPoints() const;

        //- Element connectivity
        void writeConnectivity
        (
            ensightGeoFile&,
            const word& key,
            const labelUList& idList,
            const labelUList& pointMap
        ) const;


        //- Helper: write connectivity
        void writeConnectivity
        (
            ensightGeoFile&,
            const word& key,
            const faceList&,
            const labelUList& idList,
            const labelUList& pointMap
        ) const;


public:

    //- Runtime type information
    TypeName("ensightFaces");


    // Constructors

        //- Construct part with 0-based index, description, points and faces
        //  Can skip local point renumbering when points are contiguous
        ensightPartFaces
        (
            label partIndex,
            const string& description,
            const pointField&,
            const faceList&,
            const bool contiguousPoints = false
        );

        //- Construct from polyMesh and polyPatch
        ensightPartFaces
        (
            label partIndex,
            const polyMesh&,
            const polyPatch&
        );

        //- Disallow default bitwise copy construction
        ensightPartFaces(const ensightPartFaces&) = delete;


    //- Destructor
    virtual ~ensightPartFaces();


    // Member Functions

        // Access

            //- Part index (0-based)
            virtual label index() const
            {
                return ensightFaces::index();
            }


            //- Number of elements in this part
            virtual label size() const
            {
                return ensightFaces::size();
            }


            //- Return the patch index, -1 when not in use.
            inline label patchIndex() const
            {
                return patchIndex_;
            }


        // Output

            //- Write summary information about the object
            virtual void writeSummary(Ostream&) const;

            //- Write geometry
            virtual void write(ensightGeoFile&) const;

            //- Helper: write geometry given the pointField
            virtual void write(ensightGeoFile&, const pointField&) const;


            //- Print various types of debugging information
            virtual void dumpInfo(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ensightPartFaces&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
