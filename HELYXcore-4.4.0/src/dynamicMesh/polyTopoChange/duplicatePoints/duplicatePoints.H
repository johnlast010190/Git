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
    (c) 2024 Engys Ltd.

Class
    Foam::duplicatePoints

Description
    Duplicate points

SourceFiles
    duplicatePoints.C

\*---------------------------------------------------------------------------*/

#ifndef duplicatePoints_H
#define duplicatePoints_H

#include "primitives/ints/lists/labelList.H"
#include "db/typeInfo/typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyTopoChange;
class edge;
class face;
class polyTopoChangeMap;
class localPointRegion;

/*---------------------------------------------------------------------------*\
                           Class duplicatePoints Declaration
\*---------------------------------------------------------------------------*/

class duplicatePoints
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Per regionSide.pointRegion() the points created.
        //  element 0 is the original point.
        //  Can have empty entries if a point only has one region.
        labelListList duplicates_;


public:

    //- Runtime type information
    ClassName("duplicatePoints");


    // Constructors

        //- Construct from mesh
        duplicatePoints(const polyMesh& mesh);

        //- Disallow default bitwise copy construction
        duplicatePoints(const duplicatePoints&) = delete;


    // Member Functions

        // Access

            //- Per point in regionSide.pointRegions() the points created.
            const labelListList& duplicates() const
            {
                return duplicates_;
            }


        // Topology changes

            //- Play commands into polyTopoChange to duplicate points. Gets
            //  localPointRegion structure which is per non-manifold point
            //  the regions per point.
            void setRefinement
            (
                const localPointRegion& regionSide,
                polyTopoChange&
            );

            //- Force recalculation of locally stored data on topological change
            void topoChange(const polyTopoChangeMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const duplicatePoints&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
