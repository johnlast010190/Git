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
    (c) 2017-2024 Engys Ltd.

Class
    Foam::GIBIntersectionData

Description
    Check if points are crossing processor boundaries and are hitting the
    boundaries. If a pointHit is found for a point then snap on the boundary.

SourceFiles
    GIBIntersectionData.C

\*---------------------------------------------------------------------------*/

#ifndef GIBIntersectionData_H
#define GIBIntersectionData_H

#include "GIB/intersectionData/intersectionData.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class GIBIntersectionData Declaration
\*---------------------------------------------------------------------------*/

class GIBIntersectionData
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Base mesh points
        const pointField& basePoints_;

        //- Base mesh face centers
        const vectorField& baseCf_;

        //- Problematic (crossing boundaries) point data in each processor
        const List<intersectionData>& iDl_;

        //- Reconstructed list which gets updated until all points investigated
        List<List<intersectionData>> giDl_;


    // Private Member Functions

        //- Intersection-checking algorithm.
        //  face --> cell --> face ... walk until it hits a boundary.
        void faceCellFaceIntersections
        (
            intersectionData&,
            LIFOStack<label>&,
            boolList&,
            boolList&
        );

        // Prepare the intersection checking algorithm
        void intersectionChecking();

        // Send problematic points
        void sendPointsToProcessors
        (
            DynamicList<intersectionData>&,
            DynamicList<label>&
        );

        // Update points
        void receiveAndUpdateData
        (
            const List<intersectionData>&,
            const List<label>&
        );

        // Check if there are any problematic points in every processor
        bool doLoop();


public:

    // Constructors

        //- Construct from components
        GIBIntersectionData
        (
            const fvMesh&,
            const pointField&,
            const vectorField&,
            const List<intersectionData>&
        );


    //- Destructor
    virtual ~GIBIntersectionData(){};


    // Member Functions

        //- Update the point displacement
        void modifyPointsPassingBoundary(pointField&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
