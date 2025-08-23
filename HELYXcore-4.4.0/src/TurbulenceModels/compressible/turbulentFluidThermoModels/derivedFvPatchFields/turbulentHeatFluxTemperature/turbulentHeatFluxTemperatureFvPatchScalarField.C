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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2010-2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "turbulentFluidThermoModels/derivedFvPatchFields/turbulentHeatFluxTemperature/turbulentHeatFluxTemperatureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "radiationModels/fvDOM/fvDOM/fvDOM.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    // declare specialization within 'Foam' namespace
    template<>
    const char* NamedEnum
    <
        Foam::compressible::
        turbulentHeatFluxTemperatureFvPatchScalarField::heatSourceType,
        4
    >::names[] =
    {
        "power",
        "flux",
        "powerTV",
        "fluxTV"
    };
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const NamedEnum
<
    turbulentHeatFluxTemperatureFvPatchScalarField::heatSourceType,
    4
> turbulentHeatFluxTemperatureFvPatchScalarField::heatSourceTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentHeatFluxTemperatureFvPatchScalarField::
turbulentHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    boundaryKappa(db(), patch(), "undefined", "undefined", "undefined-K"),
    heatSource_(hsPower),
    q_(p.size(), 0.0),
    qtimevar_(),
    QrName_("undefinedQr"),
    relax_(1.0)
{}


turbulentHeatFluxTemperatureFvPatchScalarField::
turbulentHeatFluxTemperatureFvPatchScalarField
(
    const turbulentHeatFluxTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    boundaryKappa(db(), patch(), ptf),
    heatSource_(ptf.heatSource_),
    q_(mapper(ptf.q_)),
    qtimevar_(ptf.qtimevar_.clone()),
    QrName_(ptf.QrName_),
    relax_(ptf.relax_)
{}


turbulentHeatFluxTemperatureFvPatchScalarField::
turbulentHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    boundaryKappa(db(), patch(), dict),
    heatSource_(heatSourceTypeNames_.read(dict.lookup("heatSource"))),
    q_(p.size(), 0.0),
    qtimevar_(),
    QrName_(dict.lookupOrDefault<word>("Qr", "none")),
    relax_(dict.lookupOrDefault<scalar>("relax", 1.0))
{

    if (dict.found("q"))
    {
        q_ = scalarField("q", dict, p.size());
    }
    else if (dict.found("qtime"))
    {
        qtimevar_ = Function1<scalar>::New("qtime", dict);
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Please supply either 'q' or"
            << "qtime" << exit(FatalIOError);
    }


    if (dict.found("gradient"))
    {
        gradient() = Field<scalar>("gradient", dict, p.size());
    }
    else
    {
        gradient() = 0.0;
    }

    // Still reading so cannot yet evaluate. Make up a value.
    fvPatchField<scalar>::operator=
    (
        patchInternalField()+gradient()*(patch().delta() & patch().nf())
    );
}


turbulentHeatFluxTemperatureFvPatchScalarField::
turbulentHeatFluxTemperatureFvPatchScalarField
(
    const turbulentHeatFluxTemperatureFvPatchScalarField& thftpsf
)
:
    fixedGradientFvPatchScalarField(thftpsf),
    boundaryKappa(db(), patch(), thftpsf),
    heatSource_(thftpsf.heatSource_),
    q_(thftpsf.q_),
    qtimevar_(thftpsf.qtimevar_.clone()),
    QrName_(thftpsf.QrName_),
    relax_(thftpsf.relax_)
{}


turbulentHeatFluxTemperatureFvPatchScalarField::
turbulentHeatFluxTemperatureFvPatchScalarField
(
    const turbulentHeatFluxTemperatureFvPatchScalarField& thftpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(thftpsf, iF),
    boundaryKappa(db(), patch(), thftpsf),
    heatSource_(thftpsf.heatSource_),
    q_(thftpsf.q_),
    qtimevar_(thftpsf.qtimevar_.clone()),
    QrName_(thftpsf.QrName_),
    relax_(thftpsf.relax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentHeatFluxTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
    m(q_, q_);
}


void turbulentHeatFluxTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const turbulentHeatFluxTemperatureFvPatchScalarField& thftptf =
        refCast<const turbulentHeatFluxTemperatureFvPatchScalarField>
        (
            ptf
        );

    q_.rmap(thftptf.q_, addr);
}


void turbulentHeatFluxTemperatureFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedGradientFvPatchScalarField::autoMapGIB(mapper);
    mapper.map(q_, scalar(0));
}


void turbulentHeatFluxTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalarField Qt(patch().size());

    switch (heatSource_)
    {
        case hsPower:
        {
            Qt = q_;
            const scalar Ap = gSum(patch().magSf());
            Qt /= Ap;
            break;
        }
        case hsFlux:
        {
            Qt = q_;
            break;
        }
        case hsTVPower:
        {
            const scalar t = db().time().timeOutputValue();
            scalar qtimevar = qtimevar_->value(t);
            Qt = qtimevar;

            const scalar Ap = gSum(patch().magSf());
            Qt /= Ap;
            break;
        }
        case hsTVFlux:
        {
            const scalar t = db().time().timeOutputValue();
            scalar qtimevar = qtimevar_->value(t);
            Qt = qtimevar;
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown heat source type. Valid types are: "
                << heatSourceTypeNames_ << nl << exit(FatalError);
        }
    }

    // Update option sources, using Newton's method to account for
    // possible face value dependence
    // We are solving
    // kappaEff*deltaCoeff*(Twall-Tpi) = Qt + boundarySource(Twall)
    // or f := Qt + boundarySource(Twall) - kappaEff*deltaCoeffs*(Twall-Tpi) = 0

    const scalarField kappaw(kappa());
    scalarField Tpi( this->patchInternalField() );
    const scalarField& deltaCoeffs = patch().deltaCoeffs();

    scalarField C1( Qt + kappaw*Tpi*deltaCoeffs );
    scalarField C2( kappaw*deltaCoeffs );

    label i = 0;
    scalar Tw_err = GREAT;
    scalar maxErr = 10*SMALL;
    label maxLoops = 100;

    scalarField Tw_new = *this;
    scalarField Tw_old;

    scalarField Qopt;
    do
    {
        Tw_old = Tw_new;
        scalarField bSourceDeriv;
        Qopt = boundarySources(Tw_old, bSourceDeriv);

        scalarField f( C1 + Qopt - C2*Tw_old );
        scalarField df( bSourceDeriv - C2 );

        Tw_new = Tw_old - f/stabilise(df, SMALL);
        Tw_err = gMax(mag((Tw_new-Tw_old)/stabilise(Tw_old, SMALL)));
        i++;
    }
    while (i < maxLoops && Tw_err > maxErr);
    // Initialise on the first iteration
    if (!QoptOld_.valid())
    {
        QoptOld_.reset(new scalarField(Qopt));
    }

    if (i == maxLoops)
    {
        WarningInFunction
            << "Non-convergence in Newton's method for "
            << "temperature-dependent boundary source, patch: "
            << this->patch().name() << nl
            << "Using old timestep for Qopt." << endl;
        Qopt = QoptOld_();
    }

    if (relax_ < 1.0)
    {
        Qopt = QoptOld_() + (1.0 - relax_)*(Qopt - QoptOld_());
    }

    QoptOld_.reset(new scalarField(Qopt));

    gradient() = (Qt+Qopt)/max(kappaw,SMALL);
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void turbulentHeatFluxTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchScalarField::write(os);
    os.writeEntry("heatSource", heatSourceTypeNames_[heatSource_]);
    boundaryKappa::write(os);
    os.writeEntry("Qr", QrName_);

    switch (heatSource_)
    {
        case hsPower:
        {
            q_.writeEntry("q", os);
            break;
        }
        case hsFlux:
        {
            q_.writeEntry("q", os);
            break;
        }
        case hsTVPower:
        {
            qtimevar_->writeData(os);
            break;
        }
        case hsTVFlux:
        {
            qtimevar_->writeData(os);
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown heat source type. Valid types are: "
                << heatSourceTypeNames_ << nl << exit(FatalError);
        }
    }

    writeEntry("value", os);

    if (relax_ < 1.0)
    {
        os.writeEntry("relax", relax_);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentHeatFluxTemperatureFvPatchScalarField
);

// Backward compatibility
addSpecialNamedToPatchFieldRunTimeSelection
(
    fvPatchScalarField,
    turbulentHeatFluxTemperatureFvPatchScalarField,
    compressible::turbulentHeatFluxTemperature,
    compressible__turbulentHeatFluxTemperature
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
