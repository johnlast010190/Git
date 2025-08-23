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
    Foam::nonConformalProcessorCyclicFvPatch

Description
    Non-conformal processor cyclic FV patch. As nonConformalCyclicFvPatch,
    but the neighbouring patch is on a different processor.

See also
    Foam::nonConformalProcessorCyclicPolyPatch
    Foam::nonConformalCyclicFvPatch

SourceFiles
    nonConformalProcessorCyclicFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalProcessorCyclicFvPatch_H
#define nonConformalProcessorCyclicFvPatch_H

#include "fvMesh/fvPatches/constraint/processorCyclic/processorCyclicFvPatch.H"
#include "fvMesh/fvPatches/constraint/nonConformalCoupled/nonConformalCoupledFvPatch.H"
#include "nonConformal/polyPatches/nonConformalProcessorCyclic/nonConformalProcessorCyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class nonConformalProcessorCyclicFvPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalProcessorCyclicFvPatch
:
    public processorCyclicFvPatch,
    public nonConformalCoupledFvPatch
{
    // Private data

        //- Reference to the polyPatch
        const nonConformalProcessorCyclicPolyPatch&
            nonConformalProcessorCyclicPolyPatch_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        virtual void makeWeights(scalarField& w) const override;


public:

    //- Runtime type information
    TypeName(nonConformalProcessorCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        nonConformalProcessorCyclicFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm
        );


    //- Destructor
    virtual ~nonConformalProcessorCyclicFvPatch();


    // Member Functions

        // Access

            //- Inherit the processor cyclic patch function
            using processorCyclicFvPatch::patch;

            //- Poly patch
            const nonConformalProcessorCyclicPolyPatch&
                nonConformalProcessorCyclicPatch() const
            {
                return nonConformalProcessorCyclicPolyPatch_;
            }

            //- Inherit the processor cyclic owner function
            using processorCyclicFvPatch::owner;

            //- Inherit the processor cyclic neighbour function
            using processorCyclicFvPatch::neighbour;

            //- Inherit the processor cyclic transform function
            using processorCyclicFvPatch::transform;

            //- Is this patch coupled? For NCC patches, the poly mesh is
            //  considered non-coupled whilst the finite volume mesh is
            //  considered coupled.
            virtual bool coupled() const override
            {
                return Pstream::parRun();
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
        virtual tmp<vectorField> delta() const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
