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
    Foam::nonConformalDiscreteMixingFvPatch

Description
    Non-conformal discrete-mixing FV patch.

See also
    Foam::nonConformalDiscreteMixingPolyPatch

SourceFiles
    nonConformalDiscreteMixingFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalDiscreteMixingFvPatch_H
#define nonConformalDiscreteMixingFvPatch_H

#include "fvMesh/fvPatches/basic/coupled/coupledFvPatch.H"
#include "fvMesh/fvPatches/constraint/nonConformalCoupled/nonConformalCoupledFvPatch.H"
#include "nonConformal/GAMG/nonConformalDiscreteMixing/nonConformalDiscreteMixingGAMGInterface.H"
#include "nonConformal/polyPatches/nonConformalDiscreteMixing/nonConformalDiscreteMixingPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class nonConformalDiscreteMixingFvPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalDiscreteMixingFvPatch
:
    public coupledFvPatch,
    public nonConformalCoupledFvPatch,
    public nonConformalLduInterface
{
    // Private data

        //- Reference to the polyPatch
        const nonConformalDiscreteMixingPolyPatch&
            nonConformalDiscreteMixingPolyPatch_;


    // Private typedefs

        //- Alias the patchToPatch intersection couple struct
        using couple = patchToPatches::intersection::couple;


    // Private member functions

        //- Make patch weighting factors
        virtual void makeWeights(scalarField&) const override;


public:

    //- Runtime type information
    TypeName(nonConformalDiscreteMixingPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch and fvBoundaryMesh
        nonConformalDiscreteMixingFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm
        );


    //- Destructor
    virtual ~nonConformalDiscreteMixingFvPatch();


    // Member Functions

        // Access

            //- Does this side own the patch ?
            virtual bool owner() const override
            {
                return nonConformalDiscreteMixingPolyPatch_.owner();
            }

            //- Poly patch
            const nonConformalDiscreteMixingPolyPatch&
            nonConformalPolyPatch() const
            {
                return nonConformalDiscreteMixingPolyPatch_;
            }

            //- Return neighbour patch ID
            virtual label nbrPatchID() const override
            {
                return nonConformalDiscreteMixingPolyPatch_.nbrPatchID();
            }

            //- Neighbour patch
            virtual const
            nonConformalDiscreteMixingFvPatch& nbrPatch() const override
            {
                return refCast<const nonConformalDiscreteMixingFvPatch>
                (
                    boundaryMesh()
                    [
                        nonConformalDiscreteMixingPolyPatch_.nbrPatchID()
                    ]
                );
            }

            //- Return the intersection engine
            virtual const
            discreteMixingIntersection& intersection() const override
            {
                return nonConformalDiscreteMixingPolyPatch_.intersection();
            }

            //- Return transformation between the coupled patches
            virtual const transformer& transform() const override
            {
                return nonConformalDiscreteMixingPolyPatch_.transform();
            }

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const override;

            //- Is this patch coupled? Returns true. For non-conformal discrete
            //  mixing patches, the poly-mesh is considered non-coupled,
            //  while the finite volume mesh is considered coupled.
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


        // Interpolations

            //- Interpolate field
            template<class Type>
            tmp<Field<Type>> interpolate(const Field<Type>& fld) const;

            //- Interpolate field with supplied operation
            template<class Type, class CombineOp>
            void interpolate
            (
                const UList<Type>& fld,
                List<Type>& result,
                const CombineOp& cop
            ) const;


       // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const override;

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData
            ) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMesh/fvPatches/constraint/nonConformalDiscreteMixing/nonConformalDiscreteMixingFvPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
