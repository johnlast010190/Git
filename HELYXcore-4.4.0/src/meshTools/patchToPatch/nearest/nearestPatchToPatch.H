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
    Foam::patchToPatches::nearest

Description
    Class to generate patchToPatch coupling geometry. Couples a face to the
    single nearest opposite face only.

SourceFiles
    nearestPatchToPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nearestPatchToPatch_H
#define nearestPatchToPatch_H

#include "patchToPatch/nearby/nearbyPatchToPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace patchToPatches
{

/*---------------------------------------------------------------------------*\
                           Class nearest Declaration
\*---------------------------------------------------------------------------*/

class nearest
:
    public nearby
{
    // Private Member Data

        //- For each source face, the distance to its coupled target face
        List<scalar> srcDistances_;

        //- For each target face, the distance to its coupled source face
        List<scalar> tgtDistances_;


    // Private Member Functions

        //- Intersect two faces
        virtual bool intersectFaces
        (
            const primitiveOldTimePatch& srcPatch,
            const vectorField& srcPointNormals,
            const vectorField& srcPointNormals0,
            const primitiveOldTimePatch& tgtPatch,
            const label srcFacei,
            const label tgtFacei
        ) override;

        //- Initialisation
        virtual void initialise
        (
            const primitiveOldTimePatch& srcPatch,
            const vectorField& srcPointNormals,
            const vectorField& srcPointNormals0,
            const primitiveOldTimePatch& tgtPatch
        ) override;

        //- Finalise the intersection locally. Trims the local target patch
        //  so that only parts that actually intersect the source remain.
        //  Returns new-to-old map for trimming target-associated data.
        virtual labelList finaliseLocal
        (
            const primitiveOldTimePatch& srcPatch,
            const vectorField& srcPointNormals,
            const vectorField& srcPointNormals0,
            const primitiveOldTimePatch& tgtPatch
        ) override;

        //- Send data that resulted from an intersection between the source
        //  patch and a distributed source-local-target patch back to the
        //  original target processes
        virtual void rDistributeTgt
        (
            const primitiveOldTimePatch& tgtPatch
        ) override;

        //- Finalise the intersection
        virtual label finalise
        (
            const primitiveOldTimePatch& srcPatch,
            const vectorField& srcPointNormals,
            const vectorField& srcPointNormals0,
            const primitiveOldTimePatch& tgtPatch,
            const transformer& tgtToSrc
        ) override;

        //- For each source face, the coupled target weights
        virtual tmpNrc<List<DynamicList<scalar>>> srcWeights() const override;

        //- For each target face, the coupled source weights
        virtual tmpNrc<List<DynamicList<scalar>>> tgtWeights() const override;


public:

    //- Runtime type information
    TypeName("nearest");


    // Constructors

        //- Construct from components
        nearest(const bool reverse);


    //- Destructor
    ~nearest();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchToPatches
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
