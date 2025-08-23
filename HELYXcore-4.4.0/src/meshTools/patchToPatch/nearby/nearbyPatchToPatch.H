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
    (c) 2022 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

Class
    Foam::patchToPatches::nearby

Description
    Class to generate patchToPatch coupling geometry. Couples a face to the
    single nearby opposite face only.

SourceFiles
    nearbyPatchToPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nearbyPatchToPatch_H
#define nearbyPatchToPatch_H

#include "patchToPatch/patchToPatch/patchToPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace patchToPatches
{

/*---------------------------------------------------------------------------*\
                            Class nearby Declaration
\*---------------------------------------------------------------------------*/

class nearby
:
    public patchToPatch
{
protected:

    // Protected Member Data

        //- For each source face, the bounding sphere
        List<Tuple2<point, scalar>> srcSpheres_;

        //- For each target face, the bounding sphere
        List<Tuple2<point, scalar>> tgtSpheres_;


    // Protected Member Functions

        //- Get the bound box for a source face
        virtual treeBoundBox srcBox
        (
            const face& srcFace,
            const pointField& srcPoints,
            const vectorField& srcPointNormals
        ) const override;

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


public:

    //- Runtime type information
    TypeName("nearby");


    // Constructors

        //- Construct from components
        nearby(const bool reverse);


    //- Destructor
    ~nearby();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchToPatches
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
