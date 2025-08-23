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
    (c) 2020 ENERCON GmbH
    (c) 2020-2022 OpenCFD Ltd.
    (c) 2024 Engys Ltd.
\*---------------------------------------------------------------------------*/

#include "derivedFvPatchFields/wallFunctions/atmEpsilonWallFunction/atmEpsilonWallFunctionFvPatchScalarField.H"
#include "derivedFvPatchFields/wallFunctions/nutWallFunctions/nutWallFunction/nutWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::atmEpsilonWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& epsilon0
)
{
    const label patchi = patch.index();
    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const scalarField& nutw = tnutw();
    const scalarField& y = turbModel.y()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradUw(mag(Uw.snGrad()));
    const scalar Cmu25 = pow025(this->Cmu());
    const scalar Cmu75 = pow(this->Cmu(), 0.75);
    const scalar kappa = this->kappa();
    const scalar yPlusLam = this->yPlusLam();
    const scalarField z0(z0_());


    #ifdef FULLDEBUG
    for (const auto& z : z0)
    {
        if (z < VSMALL)
        {
            FatalErrorInFunction
                << "z0 field can only contain positive values. "
                << "Please check input field z0."
                << exit(FatalError);
        }
    }
    #endif

    const labelUList& faceCells = patch.faceCells();

    // Set epsilon and G
    forAll(nutw, facei)
    {
        const label celli = faceCells[facei];

        const scalar yPlus = Cmu25*y[facei]*sqrt(k[celli])/nuw[facei];

        const scalar w = cornerWeights[facei];

        // (PGVB:Eq. 7, RH:Eq. 8)
        scalar epsilonc =
            w*Cmu75*pow(k[celli], 1.5)/(kappa*(y[facei] + z0[facei]));

        scalar Gc =
            w
            *(nutw[facei] + nuw[facei])
            *magGradUw[facei]
            *Cmu25*sqrt(k[celli])
            /(kappa*(y[facei] + z0[facei]));

        if (lowReCorrection_ && yPlus < yPlusLam)
        {
            epsilonc = w*2.0*k[celli]*nuw[facei]/sqr(y[facei] + z0[facei]);
            Gc = 0;
        }

        epsilon0[celli] += epsilonc;

        G0[celli] += Gc;
    }
}


void Foam::atmEpsilonWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    os.writeEntryIfDifferent<bool>("lowReCorrection", false, lowReCorrection_);

    if (z0_.valid())
    {
        z0_->writeEntry("roughnessHeight", os);
    }
    epsilonWallFunctionFvPatchScalarField::writeLocalEntries(os);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::atmEpsilonWallFunctionFvPatchScalarField::
atmEpsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    epsilonWallFunctionFvPatchScalarField(p, iF),
    z0_(nullptr)
{}


Foam::atmEpsilonWallFunctionFvPatchScalarField::
atmEpsilonWallFunctionFvPatchScalarField
(
    const atmEpsilonWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    epsilonWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    z0_()
{
    if (ptf.z0_.valid())
    {
        z0_.reset(mapper(ptf.z0_()).ptr());
    }
}


Foam::atmEpsilonWallFunctionFvPatchScalarField::
atmEpsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    epsilonWallFunctionFvPatchScalarField(p, iF, dict),
    z0_(nullptr)
{
    if (dict.found("roughnessHeight"))
    {
        z0_.reset
        (
            new scalarField("roughnessHeight", dict, p.size())
        );
    }
    else
    {
        z0_.reset
        (
            new scalarField(p.size(), 0)
        );
    }
}


Foam::atmEpsilonWallFunctionFvPatchScalarField::
atmEpsilonWallFunctionFvPatchScalarField
(
    const atmEpsilonWallFunctionFvPatchScalarField& ewfpsf
)
:
    epsilonWallFunctionFvPatchScalarField(ewfpsf),
    z0_()
{
    if (ewfpsf.z0_.valid())
    {
        z0_.reset(new scalarField(ewfpsf.z0_()));
    }
}


Foam::atmEpsilonWallFunctionFvPatchScalarField::
atmEpsilonWallFunctionFvPatchScalarField
(
    const atmEpsilonWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    epsilonWallFunctionFvPatchScalarField(ewfpsf, iF),
    z0_()
{
    if (ewfpsf.z0_.valid())
    {
        z0_.reset(new scalarField(ewfpsf.z0_()));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::atmEpsilonWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    epsilonWallFunctionFvPatchScalarField::autoMap(m);
    if (z0_.valid())
    {
        m(*z0_, *z0_);
    }
}


void Foam::atmEpsilonWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    epsilonWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const auto& atmpsf =
        refCast<const atmEpsilonWallFunctionFvPatchScalarField>(ptf);
    if (z0_.valid())
    {
        z0_->rmap(atmpsf.z0_(), addr);
    }
}

void Foam::atmEpsilonWallFunctionFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    epsilonWallFunctionFvPatchScalarField::autoMapGIB(mapper);
    if (z0_.valid())
    {
        mapper.map(z0_(), scalar(0));
    }
}


void Foam::atmEpsilonWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    writeLocalEntries(os);
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        atmEpsilonWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
