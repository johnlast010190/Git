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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2017-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "meshes/polyMesh/polyDistributionMap/distributionMap.H"
#include "derivedFvPatchFields/phaseChangeHumidity/phaseChangeHumidityFvPatchScalarField.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"
#include "turbulentFluidThermoModels/derivedFvPatchFields/thermalBaffle1D/thermalBaffle1DFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressible::thermalBaffle1DFvPatchScalarField::
thermalBaffle1DFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mappedPatchBase(p.patch()),
    mixedFvPatchScalarField(p, iF),
    TName_("T"),
    baffleActivated_(true),
    thickness_(p.size()),
    qs_(p.size()),
    solidDict_(),
    QPrevious_(p.size()),
    QRelaxation_(1)
{}


Foam::compressible::thermalBaffle1DFvPatchScalarField::
thermalBaffle1DFvPatchScalarField
(
    const thermalBaffle1DFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mappedPatchBase
    (
        p.patch(),
        ptf,
        isA<directFvPatchFieldMapper>(mapper)
      ? dynamic_cast<const directFvPatchFieldMapper&>(mapper).addressing()
      : dynamic_cast<const generalFvPatchFieldMapper&>(mapper)
       .directAddressing()
    ),
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    baffleActivated_(ptf.baffleActivated_),
    thickness_(mapper(ptf.thickness_)),
    qs_(mapper(ptf.qs_)),
    solidDict_(ptf.solidDict_),
    QPrevious_(mapper(ptf.QPrevious_)),
    QRelaxation_(ptf.QRelaxation_)
{}


Foam::compressible::thermalBaffle1DFvPatchScalarField::
thermalBaffle1DFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mappedPatchBase(p.patch(), dict),
    mixedFvPatchScalarField(p, iF, dict, false),
    TName_("T"),
    baffleActivated_(dict.lookupOrDefault<bool>("baffleActivated", true)),
    thickness_(),
    qs_(p.size(), 0),
    solidDict_(dict),
    QPrevious_(p.size(), 0.0),
    QRelaxation_(dict.lookupOrDefault<scalar>("relaxation", 1))
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("thickness"))
    {
        thickness_ = scalarField("thickness", dict, p.size());
    }

    if (dict.found("qs"))
    {
        qs_ = scalarField("qs", dict, p.size());
    }

    if (dict.found("QPrevious"))
    {
        QPrevious_ = scalarField("QPrevious", dict, p.size());
    }

    if (dict.found("refValue") && baffleActivated_)
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume zeroGradient.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 0.0;
    }

}


Foam::compressible::thermalBaffle1DFvPatchScalarField::
thermalBaffle1DFvPatchScalarField
(
    const thermalBaffle1DFvPatchScalarField& ptf
)
:
    mappedPatchBase(ptf.patch().patch(), ptf),
    mixedFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    baffleActivated_(ptf.baffleActivated_),
    thickness_(ptf.thickness_),
    qs_(ptf.qs_),
    solidDict_(ptf.solidDict_),
    QPrevious_(ptf.QPrevious_),
    QRelaxation_(ptf.QRelaxation_)
{}


Foam::compressible::thermalBaffle1DFvPatchScalarField::
thermalBaffle1DFvPatchScalarField
(
    const thermalBaffle1DFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mappedPatchBase(ptf.patch().patch(), ptf),
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    baffleActivated_(ptf.baffleActivated_),
    thickness_(ptf.thickness_),
    qs_(ptf.qs_),
    solidDict_(ptf.solidDict_),
    QPrevious_(ptf.QPrevious_),
    QRelaxation_(ptf.QRelaxation_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::compressible::thermalBaffle1DFvPatchScalarField::owner() const
{
    const label patchi = patch().index();

    const label nbrPatchi = samplePolyPatch().index();

    return (patchi < nbrPatchi);
}


const Foam::compressible::thermalBaffle1DFvPatchScalarField&
Foam::compressible::thermalBaffle1DFvPatchScalarField::neighbour() const
{
    return
        refCast<const thermalBaffle1DFvPatchScalarField>
        (
            patch().boundaryMesh()[samplePolyPatch().index()].template
                lookupPatchField<volScalarField, scalar>(TName_)
        );
}


const Foam::baseModels<Foam::scalar>&
Foam::compressible::thermalBaffle1DFvPatchScalarField::kappaModel() const
{
    if (this->owner())
    {
        const word phaseName(internalField().group());
        materialTables& tables =
            db().subRegistry("materialModels").lookupObjectRef<materialTables>
            (
                "materialTables"
            );
        const word kappaName(solidDict_.lookup<word>("kappa"));
        if (!tables.foundModel(tables.tableName(phaseName, ""), kappaName))
        {
            tables.addSpeciesAndSpeciesMixtures
            (
                kappaName + "Model",
                wordList({kappaName}),
                phaseName,
                word::null,
                word::null,
                wordList({kappaModel::typeName})
            );
        }
        return tables(kappaName, phaseName);
    }
    else
    {
        return neighbour().kappaModel();
    }
}


bool Foam::compressible::thermalBaffle1DFvPatchScalarField::isLegacy() const
{
    if (this->owner())
    {
        return solidDict_.found("transport");
    }
    else
    {
        return neighbour().isLegacy();
    }
}


Foam::scalar
Foam::compressible::thermalBaffle1DFvPatchScalarField::kappa() const
{
    if (this->owner())
    {
        return
            solidDict_.subDict("transport").lookup<scalar>("kappa");
    }
    else
    {
        return neighbour().kappa();
    }
}


Foam::tmp<Foam::scalarField>
Foam::compressible::thermalBaffle1DFvPatchScalarField::baffleThickness() const
{
    if (this->owner())
    {
        if (thickness_.size() != patch().size())
        {
            FatalErrorInFunction
            << " Field thickness has not been specified "
            << " for patch " << this->patch().name()
            << " in dictionary " <<  solidDict_
            << abort(FatalError);
        }

        return thickness_;
    }
    else
    {
        // Using GREAT as a default value can be problematic due to roundoff
        // in the weights; rather interpolate 1/thickness
        tmp<scalarField> trThickness =
            scalar(1)/stabilise(neighbour().baffleThickness(), SMALL);
        scalarField& rThickness = trThickness.ref();
        scalarField defaultValues(patch().size(), SMALL);
        this->distribute
        (
            rThickness,
            static_cast<const UList<scalar>&>(defaultValues)
        );
        return 1/trThickness;
    }
}


Foam::tmp<Foam::scalarField>
Foam::compressible::thermalBaffle1DFvPatchScalarField::qs() const
{
    if (this->owner())
    {
        return qs_;
    }
    else
    {
        tmp<scalarField> tqs(new scalarField(neighbour().qs()));
        scalarField& qs = tqs.ref();

        const scalarField defaultValues(patch().size(), 0.0);
        this->distribute
        (
            qs,
            static_cast<const UList<scalar>&>(defaultValues)
        );

        return tqs;
    }
}


void Foam::compressible::thermalBaffle1DFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mappedPatchBase::clearOut();

    mixedFvPatchScalarField::autoMap(m);

    if (this->owner())
    {
        m(thickness_, thickness_);
        m(qs_, qs_);
    }

    m(QPrevious_, QPrevious_);
}


void Foam::compressible::thermalBaffle1DFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const thermalBaffle1DFvPatchScalarField& tiptf =
        refCast<const thermalBaffle1DFvPatchScalarField>(ptf);

    if (this->owner())
    {
        thickness_.rmap(tiptf.thickness_, addr);
        qs_.rmap(tiptf.qs_, addr);
    }
    QPrevious_.rmap(tiptf.QPrevious_, addr);
}


void Foam::compressible::thermalBaffle1DFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    mixedFvPatchScalarField::reset(ptf);

    const thermalBaffle1DFvPatchScalarField& tiptf =
        refCast<const thermalBaffle1DFvPatchScalarField>(ptf);

    if (this->owner())
    {
        thickness_.reset(tiptf.thickness_);
        qs_.reset(tiptf.qs_);
    }
    QPrevious_.reset(tiptf.QPrevious_);
}


void Foam::compressible::thermalBaffle1DFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::autoMapGIB(mapper);
    mapper.map(thickness_, scalar(0));
    mapper.map(qs_, scalar(0));
    mapper.map(QPrevious_, scalar(0));
}


void Foam::compressible::thermalBaffle1DFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    const label patchi = patch().index();

    const label nbrPatchi = samplePolyPatch().index();

    if (baffleActivated_)
    {
        const fvPatch& nbrPatch = patch().boundaryMesh()[nbrPatchi];

        const compressible::turbulenceModel& turbModel =
            db().template lookupObject<compressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        // local properties
        const scalarField kappaw(turbModel.kappaEff(patchi));

        const fvPatchScalarField& Tp =
            patch().template lookupPatchField<volScalarField, scalar>(TName_);

        tmp<scalarField> Ti = patchInternalField();

        scalarField myh(patch().deltaCoeffs()*kappaw);

        // nbr properties
        const scalarField nbrKappaw(turbModel.kappaEff(nbrPatchi));

        scalarField nbrTi(neighbour().patchInternalField());
        // For NR solver
        scalarField nbrTi_orig(nbrTi);
        this->distribute
        (
            nbrTi,
            static_cast<const UList<scalar>&>(Ti)
        );

        scalarField nbrTp =
            neighbour().template
            lookupPatchField<volScalarField, scalar>(TName_);
        // for NR solver
        scalarField nbrTp_orig(nbrTp);
        this->distribute
        (
            nbrTp,
            static_cast<const UList<scalar>&>(Tp)
        );

        scalarField nbrh(nbrPatch.deltaCoeffs()*nbrKappaw);
        scalarField nbrh_orig(nbrh); //for NR solver
        this->distribute
        (
            nbrh,
            static_cast<const UList<scalar>&>(patch().deltaCoeffs()*kappaw)
        );

        // Calculate boundary solid thermo kappa
        const scalarField Tback(*this);
        scalarField& TpRef = *this;
        TpRef = 0.5*(Tp + nbrTp);

        // Legacy support
        scalarField kappas;
        if (isLegacy())
        {
            kappas = scalarField(this->size(), kappa());
        }
        else
        {
            kappas = scalarField(kappaModel().boundaryField()[patchi]);
        }
        TpRef = Tback;
        const scalarField KDeltaw(kappas/baffleThickness());

        scalarField Qfixed(0.5*qs_);

        // Update option sources, using Newton's method to account for
        // possible face value dependence
        // We are solving
        // KDeltaw*(Twall-nbrTp) + myh*(Twall+Ti) = Qt + boundarySource(Twall)
        scalarField Tpi(this->patchInternalField());

        scalarField C1(Qfixed + myh*Ti + KDeltaw*nbrTp);
        scalarField C2(myh + KDeltaw);

        label i = 0;
        scalar Tw_err = GREAT;
        scalar maxErr = 1e-5;
        label maxLoops = 100;

        scalarField Tw_new = *this;
        scalarField Tw_old;

        scalarField Qopt;
        do
        {
            Tw_old = Tw_new;
            scalarField bSourceDeriv;
            Qopt = boundarySources(Tw_old, bSourceDeriv);
            scalarField f(C1 + Qopt - C2*Tw_old);
            scalarField df(bSourceDeriv - C2);
            Tw_new = Tw_old - f/df;
            Tw_err = gMax(mag(Tw_new-Tw_old)/stabilise(Tw_old, SMALL));
            i++;
        }
        while (i < maxLoops && Tw_err > maxErr);

        if (i == maxLoops)
        {
            WarningInFunction
                << "Non-convergence in Newton's method for "
                << "temperature-dependent boundary source, patch: "
                << this->patch().name() << nl << endl;
        }

        scalarField Qtot(Qfixed + Qopt);
        Qtot = QRelaxation_*Qtot + (1-QRelaxation_)*QPrevious_;
        QPrevious_ = Qtot;

        valueFraction() = 1.0/(1.0 + myh/KDeltaw);
        refValue() = 1.0/C2*(KDeltaw*nbrTp + Qtot)/valueFraction();

        if (debug)
        {
            const scalarField qDot(kappaw*Tp.snGrad());
            scalar Q = gSum(patch().magSf()*qDot);
            Info<< patch().boundaryMesh().mesh().name() << ':'
                << patch().name() << ':'
                << this->internalField().name() << " <- "
                << nbrPatch.name() << ':'
                << this->internalField().name() << " :"
                << " heat[W]:" << Q
                << " walltemperature "
                << " min:" << gMin(*this)
                << " max:" << gMax(*this)
                << " avg:" << gAverage(*this)
                << endl;
        }
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::compressible::thermalBaffle1DFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    mappedPatchBase::write(os);

    if (this->owner())
    {
        baffleThickness()().writeEntry("thickness", os);
        qs()().writeEntry("qs", os);
        if (isLegacy())
        {
            os.writeEntry("transport", solidDict_.subDict("transport"));
        }
        else
        {
            os.writeEntry("kappa", solidDict_.lookup<word>("kappa"));
        }
    }

    QPrevious_.writeEntry("QPrevious", os);
    os.writeEntry("relaxation", QRelaxation_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
    makePatchTypeField(fvPatchScalarField, thermalBaffle1DFvPatchScalarField);

    // Backward compatibility
    addSpecialNamedToPatchFieldRunTimeSelection
    (
        fvPatchScalarField,
        thermalBaffle1DFvPatchScalarField,
        thermalBaffle1D<hConstSolidThermoPhysics>,
        thermalBaffle1D__legacy
    );
}
}

// ************************************************************************* //
