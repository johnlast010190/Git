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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2022-2025 Engys Ltd.

Class
    Foam::nonConformalCyclicFvPatch

Description
    Non-conformal cyclic FV patch. As nonConformalCoupledFvPatch, but the
    neighbouring patch is local and known and is made available by this class.

See also
    Foam::nonConformalCyclicPolyPatch
    Foam::nonConformalCoupledFvPatch

SourceFiles
    nonConformalCyclicFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalCyclicFvPatch_H
#define nonConformalCyclicFvPatch_H

#include "fvMesh/fvPatches/constraint/cyclic/cyclicFvPatch.H"
#include "fvMesh/fvPatches/constraint/nonConformalCoupled/nonConformalCoupledFvPatch.H"
#include "nonConformal/polyPatches/nonConformalCyclic/nonConformalCyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class nonConformalErrorFvPatch;

/*---------------------------------------------------------------------------*\
                  Class nonConformalCyclicFvPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalCyclicFvPatch
:
    public cyclicFvPatch,
    public nonConformalCoupledFvPatch
{
    // Private data

        //- Reference to the polyPatch
        const nonConformalCyclicPolyPatch& nonConformalCyclicPolyPatch_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        virtual void makeWeights(scalarField&) const override;


public:

    //- Runtime type information
    TypeName(nonConformalCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch and fvBoundaryMesh
        nonConformalCyclicFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm
        );


    //- Destructor
    virtual ~nonConformalCyclicFvPatch();


    // Member Functions

        // Access

            //- Inherit the cyclic patch function
            using cyclicFvPatch::patch;

            //- Poly patch
            const nonConformalCyclicPolyPatch& nonConformalCyclicPatch() const
            {
                return nonConformalCyclicPolyPatch_;
            }

            //- Inherit the cyclic owner function
            using cyclicFvPatch::owner;

            //- Inherit the cyclic neighbour function
            using cyclicFvPatch::neighbour;

            //- Inherit the cyclic transform function
            using cyclicFvPatch::transform;

            //- Neighbour patch
            virtual const nonConformalCyclicFvPatch& nbrPatch() const override;

            //- Is this patch coupled? Returns true. For NCC patches, the poly
            //  mesh is considered non-coupled whilst the finite volume mesh is
            //  considered coupled.
            virtual bool coupled() const override
            {
                return true;
            }

            //- Return the start label of this patch in the polyMesh face list.
            //  Raise an error unless the patch is empty.
            virtual label start() const override
            {
                return nonConformalFvPatch::start();
            }

            //- Return the size
            virtual label size() const override
            {
                return nonConformalFvPatch::size();
            }

            //- Return the face-cells
            virtual const labelUList& faceCells() const override
            {
                return nonConformalFvPatch::faceCells();
            }


        //- Return delta (P to N) vectors across coupled patch
        virtual tmp<vectorField> delta() const override
        {
            return coupledFvPatch::delta(nbrPatch().coupledFvPatch::delta());
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
