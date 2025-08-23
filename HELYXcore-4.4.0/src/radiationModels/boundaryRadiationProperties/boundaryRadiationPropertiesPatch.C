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
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2010-2016 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "boundaryRadiationProperties/boundaryRadiationPropertiesPatch.H"
#include "mappedPatches/mappedPolyPatch/mappedPatchBase.H"
#include "radiationModels/radiationModel/radiationModel.H"
#include "global/constants/physicoChemical/physicoChemicalConstants.H"
#include "absorptionEmissionModels/absorptionEmissionModel/absorptionEmissionModel.H"
#include "primitives/functions/Function1/Constant/Constant.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::radiationModels::boundaryRadiationPropertiesPatch::methodType,
        3
    >::names[] =
    {
        "solidRadiation",
        "lookup",
        "model"
    };

    template<>
    const char* Foam::NamedEnum
    <
        Foam::radiationModels::boundaryRadiationPropertiesPatch::modelType,
        2
    >::names[] =
    {
        "transparent",
        "opaque"
    };

}

const Foam::NamedEnum
<
    Foam::radiationModels::boundaryRadiationPropertiesPatch::methodType,
    3
> Foam::radiationModels::boundaryRadiationPropertiesPatch::methodTypeNames_;

const Foam::NamedEnum
<
    Foam::radiationModels::boundaryRadiationPropertiesPatch::modelType,
    2
> Foam::radiationModels::boundaryRadiationPropertiesPatch::modelTypeNames_;

// * * * * * * * * * * * * * * * * Private functions * * * * * * * * * * * * //

Foam::label
Foam::radiationModels::boundaryRadiationPropertiesPatch::nbrPatchIndex() const
{
    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch_);

    return (mpp.samplePolyPatch().index());
}


const Foam::fvMesh&
Foam::radiationModels::boundaryRadiationPropertiesPatch::nbrRegion() const
{
    const mappedPatchBase& mpp = refCast<const mappedPatchBase>(patch_);

    return (refCast<const fvMesh>(mpp.sampleMesh()));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::boundaryRadiationPropertiesPatch::
boundaryRadiationPropertiesPatch
(
    const polyPatch& p,
    const dictionary& dict,
    const scalar& TRef
)
:
    method_(methodTypeNames_.read(dict.lookup("mode"))),
    model_(modelTypeNames_[dict.lookupOrDefault<word>("modelType", "opaque")]),
    dict_(dict),
    absorptionEmission_(nullptr),
    transmissivity_(nullptr),
    Ta_(),
    Qa_(),
    patch_(p),
    TRef_(TRef)
{
    switch (method_)
    {
        case SOLIDRADIATION:
        {
            if (!isA<mappedPatchBase>(p))
            {
                FatalErrorInFunction
                    << "\n    patch type '" << p.type()
                    << "' not type '" << mappedPatchBase::typeName << "'"
                    << "\n    for patch " << p.name()
                    << abort(FatalIOError);
            }
        }
        break;

        case MODEL:
        {
            const fvMesh& mesh =
                refCast<const fvMesh>(p.boundaryMesh().mesh());

            absorptionEmission_.reset
            (
                absorptionEmissionModel::New(dict, mesh).ptr()
            );

            transmissivity_.reset
            (
                transmissivityModel::New(dict, mesh).ptr()
            );
        }
        break;

        case LOOKUP:
        {
            //Do nothing
        }
        break;
    }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::boundaryRadiationPropertiesPatch::
~boundaryRadiationPropertiesPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::radiationModels::absorptionEmissionModel&
Foam::radiationModels::boundaryRadiationPropertiesPatch::absorptionEmission() const
{
    return absorptionEmission_();
}


const Foam::radiationModels::transmissivityModel&
Foam::radiationModels::boundaryRadiationPropertiesPatch::transmissiveModel() const
{
    return transmissivity_();
}


Foam::tmp<Foam::scalarField>
Foam::radiationModels::boundaryRadiationPropertiesPatch::emissivity
(
    const label bandI
) const
{
    switch (method_)
    {
        case SOLIDRADIATION:
        {
            const fvMesh& nbrMesh = nbrRegion();

            const radiationModel& radiation =
                nbrMesh.lookupObject<radiationModel>
                (
                    "radiationProperties"
                );

            tmp<scalarField> te
            (
                new scalarField
                (
                    radiation.absorptionEmission().e(bandI)().boundaryField()
                    [
                        nbrPatchIndex()
                    ]
                )
            );

            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>(patch_);

            scalarField defaultValues
            (
                patch_.size(), dict_.lookupOrDefault("emissivity", scalar(0))
            );

            mpp.distribute(te.ref(), dynamic_cast<UList<scalar>&>(defaultValues));

            return te;

        }
        break;

        case LOOKUP:
        {
            tmp<scalarField> e
            (
                new scalarField
                (
                    patch_.size(),
                    dict_.lookup<scalar>("emissivity")
                )
            );

            return e;
        }

        case MODEL:
        {
            const label index = patch_.index();
            tmp<scalarField> e
            (
                new scalarField
                (
                   absorptionEmission_->e(bandI)().boundaryField()[index]
                )
            );

            return e;
        }

        default:
        {
            FatalErrorInFunction
                << "Please set 'mode' to one of "
                << methodTypeNames_
                << exit(FatalError);
        }
        break;
    }

    return scalarField(0);
}


Foam::tmp<Foam::scalarField>
Foam::radiationModels::boundaryRadiationPropertiesPatch::absorptivity
(
    const label bandI
) const
{
    switch (method_)
    {
        case SOLIDRADIATION:
        {
            const fvMesh& nbrMesh = nbrRegion();

            const radiationModel& radiation =
                nbrMesh.lookupObject<radiationModel>("radiationProperties");

            tmp<scalarField> ta
            (
                new scalarField
                (
                    radiation.absorptionEmission().a(bandI)().boundaryField()
                    [
                        nbrPatchIndex()
                    ]
                )
            );

            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>(patch_);

            scalarField defaultValues
            (
                patch_.size(),
                dict_.lookupOrDefault("absorptivity", scalar(0))
            );

            mpp.distribute
            (
                ta.ref(),
                dynamic_cast<UList<scalar>&>(defaultValues)
            );

            return ta;

        }
        break;

        case MODEL:
        {
            const label index = patch_.index();
            tmp<scalarField> a
            (
                new scalarField
                (
                   absorptionEmission_->a(bandI)().boundaryField()[index]
                )
            );
            return a;
        }

        case LOOKUP:
        {
            tmp<scalarField> a
            (
                new scalarField
                (
                    patch_.size(),
                    dict_.lookup<scalar>("absorptivity")
                )
            );

            return a;
        }

        default:
        {
            FatalErrorInFunction
                << "Unimplemented method " << method_ << endl
                << "Please set 'mode' to one of "
                << methodTypeNames_
                << exit(FatalError);
        }
        break;
    }

    return scalarField(0);
}


Foam::tmp<Foam::scalarField>
Foam::radiationModels::boundaryRadiationPropertiesPatch::transmissivity
(
    const label bandI,
    bool useSolarTransmissivity
) const
{
    switch (method_)
    {
        case SOLIDRADIATION:
        {
            const fvMesh& nbrMesh = nbrRegion();

            const radiationModel& radiation =
                nbrMesh.lookupObject<radiationModel>("radiationProperties");

            tmp<scalarField> tt
            (
                new scalarField
                (
                    radiation.transmissivity().tauEff(bandI)().boundaryField()
                    [
                        nbrPatchIndex()
                    ]
                )
            );

            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>(patch_);

            scalarField defaultValues
            (
                patch_.size(),
                (useSolarTransmissivity && dict_.found("transmissivitySolar"))
              ? dict_.lookup<scalar>("transmissivitySolar")
              : dict_.lookupOrDefault("transmissivity", scalar(0))
            );

            mpp.distribute
            (
                tt.ref(), dynamic_cast<UList<scalar>&>(defaultValues)
            );

            return tt;

        }
        break;

        case MODEL:
        {
            switch (model_)
            {
                case TRANSPARENT:
                {
                    return tmp<scalarField>::New(patch_.size(), 1.0);
                }
                break;

                case OPAQUE:
                {
                    return tmp<scalarField>::New(patch_.size(), 0.0);
                }
                break;
            }
        }
        break;

        case LOOKUP:
        {
            tmp<scalarField> tau
            (
                (useSolarTransmissivity && dict_.found("transmissivitySolar"))
              ? new scalarField
                (
                        patch_.size(),
                        dict_.lookup<scalar>("transmissivitySolar")
                )
              : new scalarField
                (
                    patch_.size(),
                    dict_.found("transmissivity")
                  ? dict_.lookup<scalar>("transmissivity")
                  : 0.0
                )
            );

            return tau;
        }
        break;

        default:
        {
            FatalErrorInFunction
                << "Unimplemented method " << method_ << endl
                << "Please set 'mode' to one of "
                << methodTypeNames_
                << exit(FatalError);
        }
        break;
    }

    return scalarField(0);
}


Foam::tmp<Foam::scalarField>
Foam::radiationModels::boundaryRadiationPropertiesPatch::reflectivity
(
    const label bandI
) const
{
    const tmp<scalarField> tt = transmissivity(bandI);
    const tmp<scalarField> ta = absorptivity(bandI);

    return (1.0 - tt - ta);
}


Foam::tmp<Foam::scalarField> Foam::radiationModels::
boundaryRadiationPropertiesPatch::
emittedRadiantFlux(const scalarField& T, const label bandI) const
{
    const tmp<scalarField> te = emissivity(bandI);
    return (te*constant::physicoChemical::sigma.value() * pow4(T + TRef_));
}


Foam::tmp<Foam::scalarField> Foam::radiationModels::
boundaryRadiationPropertiesPatch::
radiantTransmissionSource(const label bandI) const
{
    const tmp<scalarField> tt = transmissivity(bandI);

    scalarField Ta(patch_.size(), Tambient()+TRef_);
    scalarField Qa(patch_.size(), Qambient());

    //diffuse transmission from outside environment
    return
    (
        (
            tt*constant::physicoChemical::sigma.value() * pow4(Ta) + Qa
        )/constant::mathematical::pi
    );
}

Foam::scalar
Foam::radiationModels::boundaryRadiationPropertiesPatch::Tambient() const
{
    if (!Ta_.valid())
    {
        if (dict_.found("Ta"))
        {
            Ta_ = Function1<scalar>::New("Ta", dict_);
        }
        else if (dict_.found("Tinf"))
        {
            Ta_ = Function1<scalar>::New("Tinf", dict_);
        }
        else //default value of 300K
        {
            Ta_ = autoPtr<Function1<scalar>>
                (new Function1Types::Constant<scalar>("Ta", 300.0-TRef_));
        }
    }

    scalar t = patch_.boundaryMesh().mesh().time().timeOutputValue();

    return Ta_->value(t);
}


Foam::scalar
Foam::radiationModels::boundaryRadiationPropertiesPatch::Qambient() const
{

    if (!Qa_.valid())
    {
        if (dict_.found("Qa"))
        {
            Qa_ = Function1<scalar>::New("Qa", dict_);
        }
        else //default value of 0 W/m2
        {
            Qa_ =
                autoPtr<Function1<scalar>>
                (
                    new Function1Types::Constant<scalar>("Qa", 0.0)
                );
        }
    }

    scalar t = patch_.boundaryMesh().mesh().time().timeOutputValue();

    return Qa_->value(t);

}


void Foam::radiationModels::boundaryRadiationPropertiesPatch::write
(
    Ostream& os
) const
{
    os.writeEntry("mode", methodTypeNames_[method_]);

    switch (method_)
    {
        case MODEL:
        {
            word modelType
            (
                word(dict_.lookup("absorptionEmissionModel"))
            );

            os.writeEntry("absorptionEmissionModel", modelType);

            word modelCoeffs(modelType + word("Coeffs"));
            os.writeKeyword(modelCoeffs);
            dict_.subDict(modelCoeffs).write(os);

            modelType = word(dict_.lookup("transmissivityModel"));

            os.writeEntry("transmissivityModel", modelType);

            modelCoeffs = modelType + word("Coeffs");
            os.writeKeyword(modelCoeffs);
            dict_.subDict(modelCoeffs).write(os);

            modelType = word(dict_.lookup("modelType"));

            os.writeEntry("modelType", modelType);

            break;
        }

        case LOOKUP:
        {
            const scalarField emissivity("emissivity", dict_, patch_.size());
            emissivity.writeEntry("emissivity", os);

            const scalarField absorptivity
            (
                "absorptivity", dict_, patch_.size()
            );
            absorptivity.writeEntry("absorptivity", os);

            const scalarField transmissivity
            (
                "transmissivity", dict_, patch_.size()
            );
            transmissivity.writeEntry("transmissivity", os);

            if (dict_.found("transmissivitySolar"))
            {
                const scalarField transmissivitySolar
                (
                    "transmissivitySolar", dict_, patch_.size()
                );
                transmissivitySolar.writeEntry("transmissivitySolar", os);
            }

            break;
        }

        case SOLIDRADIATION:
        {
        }
    }
}


// ************************************************************************* //
