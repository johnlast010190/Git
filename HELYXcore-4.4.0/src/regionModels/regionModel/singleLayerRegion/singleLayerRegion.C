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
    (c) 2011-2020 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "regionModel/singleLayerRegion/singleLayerRegion.H"
#include "fvMesh/fvMesh.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "db/Time/Time.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
    defineTypeNameAndDebug(singleLayerRegion, 0);
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::singleLayerRegion::singleLayerRegion
(
    const fvMesh& mesh,
    const word& regionType,
    const word& modelName,
    bool readFields
)
:
    regionModel(mesh, regionType, modelName, false),
    nHat_
    (
        IOobject("nHat", time_.timeName(), regionMesh()),
        regionMesh(),
        dimensionedVector(dimless, Zero),
        zeroGradientFvPatchField<vector>::typeName
    ),
    magSf_
    (
        IOobject("magSf", time_.timeName(), regionMesh()),
        regionMesh(),
        dimensionedScalar(dimArea, 0)
    ),
    VbyA_
    (
        IOobject("d", time_.timeName(), regionMesh()),
        regionMesh(),
        dimensionedScalar(dimLength, 0),
        zeroGradientFvPatchField<vector>::typeName
    ),
    passivePatchIDs_()
{
    label nBoundaryFaces = 0;
    const polyBoundaryMesh& rbm = regionMesh().boundaryMesh();
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        const polyPatch& pp = rbm[patchi];
        const labelList& fCells = pp.faceCells();

        nBoundaryFaces += fCells.size();

        UIndirectList<vector>(nHat_, fCells) = pp.faceNormals();
        UIndirectList<scalar>(magSf_, fCells) = mag(pp.faceAreas());
    }
    nHat_.correctBoundaryConditions();

    if (nBoundaryFaces != regionMesh().nCells())
    {
        FatalErrorInFunction
            << "Number of primary region coupled boundary faces not equal to "
            << "the number of cells in the local region" << nl << nl
            << "Number of cells = " << regionMesh().nCells() << nl
            << "Boundary faces  = " << nBoundaryFaces << nl
            << abort(FatalError);
    }

    passivePatchIDs_.setSize(intCoupledPatchIDs_.size(), -1);
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        const polyPatch& ppIntCoupled = rbm[patchi];
        if (ppIntCoupled.size() > 0)
        {
            const label cellId = rbm[patchi].faceCells()[0];
            const cell& cFaces = regionMesh().cells()[cellId];
            const label faceO
            (
                cFaces.opposingFaceLabel
                (
                    ppIntCoupled.start(), regionMesh().faces()
                )
            );
            passivePatchIDs_[i] = rbm.whichPatch(faceO);
        }
    }

    Pstream::listCombineGather(passivePatchIDs_, maxEqOp<label>());
    Pstream::listCombineScatter(passivePatchIDs_);

    VbyA_.primitiveFieldRef() = regionMesh().V()/magSf_;
    VbyA_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionModels::singleLayerRegion::~singleLayerRegion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::volVectorField& Foam::regionModels::singleLayerRegion::nHat() const
{
    return nHat_;
}


const Foam::volScalarField::Internal&
Foam::regionModels::singleLayerRegion::magSf() const
{
    return magSf_;
}


const Foam::volScalarField& Foam::regionModels::singleLayerRegion::VbyA() const
{
    return VbyA_;
}


const Foam::labelList&
Foam::regionModels::singleLayerRegion::passivePatchIDs() const
{
    return passivePatchIDs_;
}


// ************************************************************************* //
