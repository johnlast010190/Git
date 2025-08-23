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

\*---------------------------------------------------------------------------*/

#include "fvMesh/fvPatches/constraint/nonConformalDiscreteMixing/nonConformalDiscreteMixingFvPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalDiscreteMixingFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, nonConformalDiscreteMixingFvPatch, polyPatch);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::nonConformalDiscreteMixingFvPatch::makeWeights(scalarField& w) const
{
    const nonConformalDiscreteMixingFvPatch& nbrPatch = this->nbrPatch();

    const scalarField nfDelta(coupledFvPatch::nf() & coupledFvPatch::delta());
    const scalarField nbrNfDelta
    (
        nbrPatch.nf() & nbrPatch.coupledFvPatch::delta()
    );

    tmp<scalarField> tInterpNbrNfDelta = interpolate(nbrNfDelta);
    const scalarField& interpNbrNfDelta = tInterpNbrNfDelta();

    forAll(nfDelta, facei)
    {
        scalar di = nfDelta[facei];
        scalar dni = interpNbrNfDelta[facei];

        w[facei] = dni/(di + dni + VSMALL);
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalDiscreteMixingFvPatch::nonConformalDiscreteMixingFvPatch
(
    const polyPatch& patch,
    const fvBoundaryMesh& bm
)
:
    coupledFvPatch(patch, bm),
    nonConformalCoupledFvPatch(static_cast<const fvPatch&>(*this)),
    nonConformalDiscreteMixingPolyPatch_
    (
        refCast<const nonConformalDiscreteMixingPolyPatch>(patch)
    )
{
    const fvMesh& mesh = bm.mesh();

    IOobject dictHeader
    (
        "dynamicMeshDict",
        mesh.time().constant(),
        mesh.dbDir(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (dictHeader.typeHeaderOk<IOdictionary>())
    {
        IOdictionary dict(dictHeader);

        if (dict.found("mover"))
        {
            FatalErrorInFunction
                << "Non-conformal discrete mixing patch types are not "
                << "compatible with moving meshes."
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalDiscreteMixingFvPatch::~nonConformalDiscreteMixingFvPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField>
Foam::nonConformalDiscreteMixingFvPatch::delta() const
{
    const nonConformalDiscreteMixingFvPatch& nbrPatch = this->nbrPatch();

    const vectorField pDelta(coupledFvPatch::delta());

    tmp<vectorField> tnbrPDelta =
        interpolate(nbrPatch.coupledFvPatch::delta()());
    const vectorField& nbrPDelta = tnbrPDelta();

    tmp<vectorField> tpdv(new vectorField(pDelta.size()));
    vectorField& pdv = tpdv.ref();

    forAll(pDelta, facei)
    {
        const vector& ddi = pDelta[facei];
        const vector& dni = nbrPDelta[facei];

        pdv[facei] = ddi - dni;
    }

    return tpdv;
}


Foam::tmp<Foam::labelField>
Foam::nonConformalDiscreteMixingFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField>
Foam::nonConformalDiscreteMixingFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    return nbrPatch().patchInternalField(iF);
}


// ************************************************************************* //
