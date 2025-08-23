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
    (c) 2021-2025 OpenFOAM Foundation
    (c) 2022-2025 Engys Ltd.

Class
    Foam::patchToPatch

Description
    Class to generate coupling geometry between two primitive patches

SourceFiles
    patchToPatch.C
    patchToPatchParallelOps.C

\*---------------------------------------------------------------------------*/

#ifndef patchToPatch_H
#define patchToPatch_H

#include "primitives/remote/remote.H"
#include "primitives/RefPair/RefPair.H"
#include "meshes/polyMesh/polyDistributionMap/distributionMap.H"
#include "PrimitiveOldTimePatch/primitiveOldTimePatch.H"
#include "meshes/treeBoundBox/treeBoundBox.H"
#include "patchToPatchTools/patchToPatchTools.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class patchToPatch Declaration
\*---------------------------------------------------------------------------*/

class patchToPatch
{
protected:

    // Protected Data

        //- Flag to indicate that the two patches are co-directional and
        //  that the orientation of one should therefore be reversed
        const bool reverse_;

        //- Index of the processor holding all faces of the patchToPatch, or -1
        //  if spread across multiple processors
        label singleProcess_;

        //- For each source face, the coupled local target faces
        List<DynamicList<label>> srcLocalTgtFaces_;

        //- For each target face, the coupled local source faces
        List<DynamicList<label>> tgtLocalSrcFaces_;

        //- Map from source patch faces to target-local source patch faces
        autoPtr<distributionMap> srcMapPtr_;

        //- Map from target patch faces to source-local target patch faces
        autoPtr<distributionMap> tgtMapPtr_;

        //- When running in parallel, a map from local source face index to
        //  source processor and face index
        autoPtr<List<remote>> localSrcProcFacesPtr_;

        //- When running in parallel, a map from local target face index to
        //  target processor and face index
        autoPtr<List<remote>> localTgtProcFacesPtr_;


    // Protected Member Functions

        // Searching

            //- Get the bound box for a source face
            virtual treeBoundBox srcBox
            (
                const face& srcFace,
                const pointField& srcPoints,
                const vectorField& srcPointNormals
            ) const = 0;

            //- Get the bound box for a source face
            treeBoundBox srcBox
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const label srcFacei
            ) const;

            //- Get the bound box for the source patch
            treeBoundBox srcBox
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0
            ) const;

            //- Get the bound box for a target face
            treeBoundBox tgtBox
            (
                const primitiveOldTimePatch& tgtPatch,
                const label tgtFacei
            ) const;

            //- Get the bound box for the target patch
            treeBoundBox tgtBox
            (
                const primitiveOldTimePatch& tgtPatch
            ) const;


        // Intersection

            //- Intersect two faces
            virtual bool intersectFaces
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const primitiveOldTimePatch& tgtPatch,
                const label srcFacei,
                const label tgtFacei
            ) = 0;

            //- Check whether two faces have already been intersected.
            //  If not, try to intersect them.
            bool findOrIntersectFaces
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const primitiveOldTimePatch& tgtPatch,
                const label srcFacei,
                const label tgtFacei
            );

            //- Intersect a queue of source-target face pairs.
            //  Update completion lists and form a new queue
            //  of target-source face pairs.
            label intersectPatchQueue
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const primitiveOldTimePatch& tgtPatch,
                const bool isSrc,
                const DynamicList<labelPair>& queue,
                labelList& faceComplete,
                DynamicList<labelPair>& otherQueue,
                const labelList& otherFaceComplete,
                boolList& otherFaceQueued,
                boolList& otherFaceVisited
            );

            //- Intersect the patches
            void intersectPatches
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const primitiveOldTimePatch& tgtPatch
            );


        // Parallel functionality

            //- Determine which target faces need to be sent to the source.
            //  This is done before intersection. Bound boxes are used to
            //  estimate what faces will intersect.
            labelListList tgtPatchSendFaces
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const primitiveOldTimePatch& tgtPatch
            ) const;

            //- Distribute a patch given its distribution map
            static List<remote> distributePatch
            (
                const distributionMap& map,
                const primitiveOldTimePatch& patch,
                autoPtr<PrimitiveOldTimePatch<faceList, pointField>>&
                    localPatchPtr
            );


        // Hooks

            //- Initialisation
            virtual void initialise
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const primitiveOldTimePatch& tgtPatch
            );

            //- Finalise the intersection locally. Trims the local target patch
            //  so that only parts that actually intersect the source remain.
            //  Returns new-to-old map for trimming target-associated data.
            virtual labelList finaliseLocal
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const primitiveOldTimePatch& tgtPatch
            );

            //- Distribute the source patch so that any information needed by
            //  the target is locally available
            virtual void distributeSrc(const primitiveOldTimePatch& srcPatch);

            //- Send data that resulted from an intersection between the source
            //  patch and a distributed source-local-target patch back to the
            //  original target processes
            virtual void rDistributeTgt(const primitiveOldTimePatch& tgtPatch);

            //- Finalising
            virtual label finalise
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const primitiveOldTimePatch& tgtPatch,
                const transformer& tgtToSrc
            );


        // Interpolation

            //- For each source face, the coupled target weights
            virtual tmpNrc<List<DynamicList<scalar>>> srcWeights() const = 0;

            //- For each target face, the coupled source weights
            virtual tmpNrc<List<DynamicList<scalar>>> tgtWeights() const = 0;


public:

    //- Runtime type information
    TypeName("patchToPatch");


    // Static Data Members

        //- Lower numeric limits for single and double precision
        static scalar small, vSmall, rootSmall, rootVSmall;
        static floatScalar floatSmall, floatVSmall;
        static floatScalar floatRootSmall, floatRootVSmall;
        static doubleScalar doubleSmall, doubleVSmall;
        static doubleScalar doubleRootSmall, doubleRootVSmall;


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        patchToPatch,
        bool,
        (const bool reverse),
        (reverse)
    );


    // Constructors

        //- Construct from components
        patchToPatch(const bool reverse);

        //- Disallow default bitwise copy construction
        patchToPatch(const patchToPatch&) = delete;


    //- Destructor
    virtual ~patchToPatch();


    // Selector

        //- Select from name
        static autoPtr<patchToPatch> New
        (
            const word& patchToPatchType,
            const bool reverse
        );


    // Member Functions

        // Access

            //- Flag to indicate that the two patches are co-directional and
            //  that the orientation of one should therefore be reversed
            bool reverse() const
            {
                return reverse_;
            }

            //- Index of the processor holding all faces of the patchToPatch,
            //  or -1 if spread across multiple processors
            label singleProcess() const
            {
                return singleProcess_;
            }

            //- Is this intersection on a single process?
            bool isSingleProcess() const
            {
                return singleProcess_ != -1;
            }

            //- Map from source patch faces to target-local source patch faces
            const distributionMap& srcMap() const
            {
                return srcMapPtr_();
            }

            //- Map from target patch faces to source-local target patch faces
            const distributionMap& tgtMap() const
            {
                return tgtMapPtr_();
            }

            //- For each source face, the coupled target procs and faces
            List<List<remote>> srcTgtProcFaces() const;

            //- For each target face, the coupled source procs and faces
            List<List<remote>> tgtSrcProcFaces() const;


        // Manipulation

            //- Update addressing and weights for the given patches
            void update
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const primitiveOldTimePatch& tgtPatch,
                const transformer& tgtToSrc = NullObjectRef<transformer>(),
                const ConstRefPair<string>& srcTgtPatchNames =
                    NullObjectRef<ConstRefPair<string>>(),
                const string& tgtToSrcName = NullObjectRef<string>()
            );

            //- Update addressing and weights for the given patches
            void update
            (
                const primitivePatch& srcPatch,
                const vectorField& srcPointNormals,
                const primitivePatch& tgtPatch,
                const transformer& tgtToSrc = NullObjectRef<transformer>(),
                const ConstRefPair<string>& srcTgtPatchNames =
                    NullObjectRef<ConstRefPair<string>>(),
                const string& tgtToSrcName = NullObjectRef<string>()
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const patchToPatch&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
