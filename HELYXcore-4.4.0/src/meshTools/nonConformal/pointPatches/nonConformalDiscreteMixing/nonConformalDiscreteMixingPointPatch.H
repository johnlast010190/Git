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
    (c) 2024 Engys Ltd.

Class
    Foam::nonConformalDiscreteMixingPointPatch

Description
    Non-conformal discrete mixing point patch for post-processing.

SourceFiles
    nonConformalDiscreteMixingPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalDiscreteMixingPointPatch_H
#define nonConformalDiscreteMixingPointPatch_H

#include "meshes/pointMesh/pointPatches/derived/coupled/coupledFacePointPatch.H"
#include "nonConformal/polyPatches/nonConformalDiscreteMixing/nonConformalDiscreteMixingPolyPatch.H"
#include "meshes/pointMesh/pointBoundaryMesh/pointBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class nonConformalDiscreteMixingPointPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalDiscreteMixingPointPatch
:
    public coupledFacePointPatch
{
    // Private Data

        //- Local reference cast into the non-conformal discrete mixing patch
        const nonConformalDiscreteMixingPolyPatch& nonConformalPolyPatch_;


    // Private Member Functions

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&) override
        {}

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&) override
        {}

        //- Initialise the patches for moving points
        virtual void initMovePoints
        (
            PstreamBuffers&, const pointField&
        ) override
        {}

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&) override
        {}

        //- Initialise the update of the patch topology
        virtual void initTopoChange(PstreamBuffers&) override;

        //- Update of the patch topology
        virtual void topoChange(PstreamBuffers&) override;


public:

    //- Runtime type information
    TypeName(nonConformalDiscreteMixingPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        nonConformalDiscreteMixingPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );

        //- Disallow default construct as copy
        nonConformalDiscreteMixingPointPatch
        (
            const nonConformalDiscreteMixingPointPatch&
        ) = delete;


    // Member Functions

        //- Return the underlying poly patch
        const nonConformalDiscreteMixingPolyPatch& nonConformalPolyPatch() const
        {
            return nonConformalPolyPatch_;
        }

        //- Return neighbour point patch
        const nonConformalDiscreteMixingPointPatch& nbrPatch() const
        {
            label patchi = nonConformalPolyPatch_.nbrPatchID();
            const pointPatch& pp = this->boundaryMesh()[patchi];

            return refCast<const nonConformalDiscreteMixingPointPatch>(pp);
        }


    // Member Operators

        //- Disallow default assignment
        void operator=(const nonConformalDiscreteMixingPointPatch&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
