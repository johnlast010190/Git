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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "motionDiffusivity/inverseDistance/inverseDistanceDiffusivity.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"
#include "fvMesh/wallDist/patchDistMethods/meshWave/meshWavePatchDistMethod.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(inverseDistanceDiffusivity, 0);

    addToRunTimeSelectionTable
    (
        motionDiffusivity,
        inverseDistanceDiffusivity,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inverseDistanceDiffusivity::inverseDistanceDiffusivity
(
    const fvMesh& mesh,
    Istream& mdData
)
:
    uniformDiffusivity(mesh, mdData),
    patchNames_(mdData),
    Lbuff_(0.0)
{
    if (!mdData.eof())
    {
        mdData >> Lbuff_;
    }

    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::inverseDistanceDiffusivity::~inverseDistanceDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::inverseDistanceDiffusivity::correct()
{
    const volScalarField& dist =
        wallDist::New
        (
            mesh(),
            patchDistMethods::meshWave::typeName,
            mesh().boundaryMesh().patchSet(patchNames_)
        ).y();

    // correct faceDiffusivity_ in presence of buffer
    if (Lbuff_ > VSMALL)
    {
        // get minimum distance (boundary layer)
        const scalar minDist = gMin(dist.primitiveField());

        // compute distance from buffer boundary (=offset)
        tmp<surfaceScalarField> tBuffDist =
            fvc::interpolate
            (
                dist - dimensionedScalar(dimLength, Lbuff_)
            );
        faceDiffusivity_ =
            dimensionedScalar(dimLength, 1)
          / max
            (
                tBuffDist,
                dimensionedScalar(dimLength, minDist)
            );
    }
    else
    {
        faceDiffusivity_ =
            dimensionedScalar(dimLength, 1)/fvc::interpolate(dist);
    }
}


// ************************************************************************* //
