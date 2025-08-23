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
    Foam::nonConformalDiscreteMixingPolyPatch

Description
    Non-conformal patch for discrete mixing coupling. Derived from
    nonConformalCoupledPolyPatch, with all owner and neighbour patch
    rotations accounted for in the intersection engine.

See also
    Foam::discreteMixingPlanePolyPatch

SourceFiles
    nonConformalDiscreteMixingPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalDiscreteMixingPolyPatch_H
#define nonConformalDiscreteMixingPolyPatch_H

#include "nonConformal/polyPatches/nonConformalCoupled/nonConformalCoupledPolyPatch.H"
#include "nonConformal/polyPatches/nonConformalDiscreteMixing/intersection/discreteMixingIntersection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class nonConformalDiscreteMixingPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalDiscreteMixingPolyPatch
:
    public coupledPolyPatch,
    public nonConformalCoupledPolyPatch
{
    // Private data

        //- Name of neighbour patch
        mutable word nbrPatchName_;

        //- Index of neighbour patch
        mutable label nbrPatchID_;

        //- Axis of rotation
        vector rotationAxis_;

        //- Centre of rotation
        point rotationCentre_;

        //- Rotation angles for the patch pair
        Pair<scalar> sectorAngles_;

        //- Number of rotational sectors for this patch
        label nSectors_;

        //- Maximum accepted absolute tolerance for the non-integer part of
        //  the number of sectors
        static constexpr scalar tolIntSec_ = 0.1;

        //- List of transforms for this patch
        List<transformer> patchTransforms_;

        //- Is the intersection engine up to date?
        mutable bool intersectionsAreValid_;

        //- Patch-to-patch intersection engine
        mutable autoPtr<discreteMixingIntersection> intersectionsPtr_;


    // Private Member functions

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&) override;

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&) override;

        //- Initialise the patches for moving points
        virtual void initMovePoints
        (
            PstreamBuffers&,
            const pointField&
        ) override
        {
            // This patch type is not compatible with moving meshes
            NotImplemented;
        }

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&) override
        {
            // This patch type is not compatible with moving meshes
            NotImplemented;
        }

        //- Initialise the update of the patch topology
        virtual void initTopoChange(PstreamBuffers&) override;

        //- Update of the patch topology
        virtual void topoChange(PstreamBuffers&) override;

        //- Reset the patch name
        virtual void rename(const wordList& newNames) override;

        //- Reset the patch index
        virtual void reorder(const labelUList& newToOldIndex) override;

        //- Compute the number of rotational sectors
        void calcSectors();

        //- Set all rotational transforms for this patch
        void setPatchTransforms();


public:

    //- Runtime type information
    TypeName("nonConformalDiscreteMixing");


    // Constructors

        //- Construct from components
        nonConformalDiscreteMixingPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from dictionary
        nonConformalDiscreteMixingPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from discrete-mixing patch components
        nonConformalDiscreteMixingPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const word& nbrPatchName,
            const word& origPatchName,
            const dictionary& discMixingDict
        );

        //- Construct as copy, resetting the boundary mesh.
        nonConformalDiscreteMixingPolyPatch
        (
            const nonConformalDiscreteMixingPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        nonConformalDiscreteMixingPolyPatch
        (
            const nonConformalDiscreteMixingPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& nbrPatchName,
            const word& origPatchName
        );

        //- Construct and return a clone, resetting the boundary mesh.
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm
        ) const override
        {
            return autoPtr<polyPatch>
            (
                new nonConformalDiscreteMixingPolyPatch(*this, bm)
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh.
        virtual autoPtr<directPolyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const override
        {
            return autoPtr<directPolyPatch>
            (
                new nonConformalDiscreteMixingPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    nbrPatchName(),
                    origPatchName()
                )
            );
        }


    //- Destructor
    virtual ~nonConformalDiscreteMixingPolyPatch();


    // Member Functions

        //- Neighbour patch name
        const word& nbrPatchName() const
        {
            return nbrPatchName_;
        }

        //- Neighbour patchID
        label nbrPatchID() const;

        //- Does this side own the patch ?
        virtual bool owner() const override
        {
            return index() < nbrPatchID();
        }

        //- Inherit the coupled neighbour method
        using coupledPolyPatch::neighbour;

        const nonConformalDiscreteMixingPolyPatch& nbrPatch() const
        {
            const polyPatch& pp = this->boundaryMesh()[nbrPatchID()];
            return refCast<const nonConformalDiscreteMixingPolyPatch>(pp);
        }

        //- Return (dummy) transformation between the coupled patches
        virtual const transformer& transform() const override
        {
            return transformer::I;
        }

        //- Return list of transforms for this patch
        const List<transformer>& patchTransforms() const
        {
            return patchTransforms_;
        }

        //- Return the intersection engine
        const discreteMixingIntersection& intersection() const
        {
            if (!owner())
            {
                FatalErrorInFunction
                    << "The non-conformal discrete mixing intersection engine "
                    << "is only available to the owner patch"
                    << abort(FatalError);
            }

            return intersectionsPtr_();
        }

        //- Access the list of patch-to-patch intersections
        const PtrList<patchToPatches::intersection>& intersections() const;

        //- Initialise ordering for primitivePatch
        virtual void initOrder
        (
            PstreamBuffers&,
            const primitivePatch&
        ) const override
        {}

        //- Return new ordering for primitivePatch
        virtual bool order
        (
            PstreamBuffers&,
            const primitivePatch&,
            labelList& faceMap,
            labelList& rotation
        ) const override
        {
            return false;
        }

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
