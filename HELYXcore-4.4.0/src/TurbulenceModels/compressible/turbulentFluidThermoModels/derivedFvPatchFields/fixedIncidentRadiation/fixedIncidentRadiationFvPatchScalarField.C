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
    (c) 2016 OpenCFD Ltd.
    (c) 2022-2024 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "turbulentFluidThermoModels/derivedFvPatchFields/fixedIncidentRadiation/fixedIncidentRadiationFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "global/constants/constants.H"
#include "radiationModels/radiationModel/radiationModel.H"
#include "absorptionEmissionModels/absorptionEmissionModel/absorptionEmissionModel.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::fixedIncidentRadiationFvPatchScalarField::
fixedIncidentRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    boundaryKappa(db(), patch(), "undefined", "undefined", "undefined-K"),
    qrIncident_(p.size(), 0.0)
{}


Foam::radiationModels::fixedIncidentRadiationFvPatchScalarField::
fixedIncidentRadiationFvPatchScalarField
(
    const fixedIncidentRadiationFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(psf, p, iF, mapper),
    boundaryKappa(db(), patch(), psf),
    qrIncident_(psf.qrIncident_)
{}


Foam::radiationModels::fixedIncidentRadiationFvPatchScalarField::
fixedIncidentRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    boundaryKappa(db(), patch(), dict),
    qrIncident_("qrIncident", dict, p.size())
{
    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=(Field<scalar>("value", dict, p.size()));
        gradient() = Field<scalar>("gradient", dict, p.size());
    }
    else
    {
        // Still reading so cannot yet evaluate. Make up a value.
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::radiationModels::fixedIncidentRadiationFvPatchScalarField::
fixedIncidentRadiationFvPatchScalarField
(
    const fixedIncidentRadiationFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(psf, iF),
    boundaryKappa(db(), patch(), psf),
    qrIncident_(psf.qrIncident_)
{}


Foam::radiationModels::fixedIncidentRadiationFvPatchScalarField::
fixedIncidentRadiationFvPatchScalarField
(
    const fixedIncidentRadiationFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    boundaryKappa(db(), patch(), ptf),
    qrIncident_(ptf.qrIncident_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiationModels::fixedIncidentRadiationFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
    m(qrIncident_, qrIncident_);
}


void Foam::radiationModels::fixedIncidentRadiationFvPatchScalarField::rmap
(
    const fvPatchScalarField& psf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(psf, addr);

    const fixedIncidentRadiationFvPatchScalarField& thftpsf =
        refCast<const fixedIncidentRadiationFvPatchScalarField>
        (
            psf
        );

    qrIncident_.rmap(thftpsf.qrIncident_, addr);
}


void
Foam::radiationModels::fixedIncidentRadiationFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedGradientFvPatchScalarField::autoMapGIB(mapper);
    mapper.map(qrIncident_, scalar(0));
}


void
Foam::radiationModels::fixedIncidentRadiationFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalarField intFld(patchInternalField());

    const radiationModel& radiation =
        db().lookupObject<radiationModel>("radiationProperties");

    scalarField emissivity
    (
        radiation.absorptionEmission().e()().boundaryField()[patch().index()]
    );

    gradient() =
        emissivity
       *(
            qrIncident_
          - physicoChemical::sigma.value()*pow4(*this)
        )/kappa();

    fixedGradientFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar qr = gSum(kappa()*gradient()*patch().magSf());
        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " -> "
            << " radiativeFlux:" << qr
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void Foam::radiationModels::fixedIncidentRadiationFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchScalarField::write(os);
    boundaryKappa::write(os);
    qrIncident_.writeEntry("qrIncident", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedIncidentRadiationFvPatchScalarField
    );
}
}

// ************************************************************************* //
