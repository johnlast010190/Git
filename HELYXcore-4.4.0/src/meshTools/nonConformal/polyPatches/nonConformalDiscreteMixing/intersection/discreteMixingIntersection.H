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
    (c) 2024-2025 Engys Ltd.

Class
    Foam::discreteMixingIntersection

Description
    Engine to generate patch-to-patch geometric intersection between two
    non-conformal discrete mixing patches.

SourceFiles
    discreteMixingIntersection.C

\*---------------------------------------------------------------------------*/

#ifndef discreteMixingIntersection_H
#define discreteMixingIntersection_H

#include "patchToPatch/intersection/intersectionPatchToPatch.H"
#include "primitives/remote/remote.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class nonConformalDiscreteMixingPolyPatch;

/*---------------------------------------------------------------------------*\
                  Class discreteMixingIntersection Declaration
\*---------------------------------------------------------------------------*/

class discreteMixingIntersection
{
    // Private Data

        //- List of patch-to-patch intersections
        PtrList<patchToPatches::intersection> patchToPatchList_;

        //- List of transforms for all combinations between owner and neighbour
        List<Pair<transformer>> transforms_;

        //- Source face areas
        scalarField srcAreas_;

        //- Target face areas
        scalarField tgtAreas_;

        //- Coupled target processor and faces across all rotational
        //  intersections, for each owner patch face.
        List<List<remote>> srcAddress_;

        //- Coupled source processor and faces across all rotational,
        //  intersections for each neighbour patch face.
        List<List<remote>> tgtAddress_;

        //- Coupled target weights across all rotational intersections,
        //  for each owner patch face.
        List<scalarList> srcWeights_;

        //- Coupled source weights across all rotational intersections,
        //  for each neighbour patch face.
        List<scalarList> tgtWeights_;

        //- Source faces with no intersections across all rotational couplings.
        //  A stabilisation procedure will be applied for these faces.
        labelHashSet srcZeroWeightFaces_;

        //- Target faces with no intersections across all rotational couplings.
        //  A stabilisation procedure will be applied for these faces.
        labelHashSet tgtZeroWeightFaces_;

        //- The intersection index of the coupled target faces,
        //  for each original owner patch face.
        List<labelList> srcIntersectionMap_;

        //- The intersection index of the coupled source faces,
        //  for each original neighbour patch face.
        List<labelList> tgtIntersectionMap_;


    // Private Member Functions

        //- Set addressing across all rotational intersections
        void setAddressing();

        //- Calculate the coupled weights for each patch face
        void calcWeights(const bool isSrc);

        //- Set source and target intersection maps
        void setIntersectionMaps();

        //- Agglomerate areas, addressing and weights. Used in GAMG solver.
        void agglomerate
        (
            const scalarField& fineAreas,
            const List<List<remote>>& fineAddress,
            const List<scalarList>& fineWeights,
            const labelList& restrictAddress,
            const labelList& otherRestrictAddress,
            const bool isSrc
        );


public:

    //- Runtime type information
    TypeName("discreteMixingIntersection");


    // Constructors

        //- Construct from poly-patches
        discreteMixingIntersection
        (
            const nonConformalDiscreteMixingPolyPatch& ownerPatch,
            const nonConformalDiscreteMixingPolyPatch& neighbourPatch
        );

        //- Construct from discrete mixing intersection and face
        //  restrict addressing.
        //  Used in agglomeration of patch-to-patch intersections.
        //  Agglomeration is passed as the addressing from fine to coarse.
        discreteMixingIntersection
        (
            const discreteMixingIntersection& fineIntersections,
            const labelList& ownRestrictAddress,
            const labelList& nbrRestrictAddress
        );

        //- Disallow default bitwise copy construction
        discreteMixingIntersection(const discreteMixingIntersection&) = delete;


    //- Destructor
    virtual ~discreteMixingIntersection();


    // Member Functions

        // Access

            //- Return the list of patch-to-patch intersection engines
            const
            PtrList<patchToPatches::intersection>& patchToPatchList() const
            {
                return patchToPatchList_;
            }

            //- Return source face areas
            const scalarField& srcAreas() const
            {
                return srcAreas_;
            }

            //- Return target face areas
            const scalarField& tgtAreas() const
            {
                return tgtAreas_;
            }

            //- Return the coupled target processor and faces for each owner
            //  face, across all rotational intersections.
            const List<List<remote>>& srcAddress() const
            {
                return srcAddress_;
            }

            //- Return the coupled source processor and faces for each neighbour
            //  face, across all rotational intersections.
            const List<List<remote>>& tgtAddress() const
            {
                return tgtAddress_;
            }

            //- Return the coupled target weights for each owner face,
            //  across all rotational intersections.
            const List<scalarList>& srcWeights() const
            {
                return srcWeights_;
            }

            //- Return the coupled source weights for each neighbour face,
            //  across all rotational intersections.
            const List<scalarList>& tgtWeights() const
            {
                return tgtWeights_;
            }

            //- Return the owner patch zero-weight faces
            const labelHashSet& srcZeroWeightFaces() const
            {
                return srcZeroWeightFaces_;
            }

            //- Return the neighbour patch zero-weight faces
            const labelHashSet& tgtZeroWeightFaces() const
            {
                return tgtZeroWeightFaces_;
            }


        // Interpolation

            //- Interpolate an owner patch field to the neighbour
            template<class Type>
            tmp<Field<Type>> interpolateToNeighbour
            (
                const Field<Type>& ownFld
            ) const;

            //- Interpolate a neighbour patch field to the owner
            template<class Type>
            tmp<Field<Type>> interpolateToOwner
            (
                const Field<Type>& nbrFld
            ) const;

            //- Interpolate to the neighbour patch with supplied operation
            template<class Type, class CombineOp>
            void interpolateToNeighbour
            (
                const UList<Type>& ownFld,
                List<Type>& result,
                const CombineOp& cop
            ) const;

            //- Interpolate to the owner patch with supplied operation
            template<class Type, class CombineOp>
            void interpolateToOwner
            (
                const UList<Type>& nbrFld,
                List<Type>& result,
                const CombineOp& cop
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const discreteMixingIntersection&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonConformal/polyPatches/nonConformalDiscreteMixing/intersection/discreteMixingIntersectionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
