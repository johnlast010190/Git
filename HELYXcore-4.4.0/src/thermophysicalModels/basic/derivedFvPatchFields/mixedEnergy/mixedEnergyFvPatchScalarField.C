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
    (c) 2011-2012 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "derivedFvPatchFields/mixedEnergy/mixedEnergyFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "basicThermo/basicThermo.H"
#include "fields/fvPatchFields/derived/mixedEnergyCalculatedTemperature/mixedEnergyCalculatedTemperatureFvPatchScalarField.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixedEnergyFvPatchScalarField::
mixedEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF)
{
    valueFraction() = 0.0;
    refValue() = 0.0;
    refGrad() = 0.0;
}


Foam::mixedEnergyFvPatchScalarField::
mixedEnergyFvPatchScalarField
(
    const mixedEnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::mixedEnergyFvPatchScalarField::
mixedEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict)
{}


Foam::mixedEnergyFvPatchScalarField::
mixedEnergyFvPatchScalarField
(
    const mixedEnergyFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf)
{}


Foam::mixedEnergyFvPatchScalarField::
mixedEnergyFvPatchScalarField
(
    const mixedEnergyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mixedEnergyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const basicThermo& thermo = basicThermo::lookupThermo(*this);
    const label patchi = patch().index();
    fvPatchScalarField& Tp =
        const_cast<fvPatchScalarField&>(thermo.T().boundaryField()[patchi]);
    const scalarField& pw = thermo.p().boundaryField()[patchi];

    if (isA<mixedFvPatchScalarField>(Tp))
    {
        mixedFvPatchScalarField& Tw = refCast<mixedFvPatchScalarField>(Tp);
        Tw.evaluate();

        valueFraction() = Tw.valueFraction();
        refValue() = thermo.he(Tw.refValue(), patchi);
        refGrad() =
            thermo.Cpv(Tw, patchi)*Tw.refGrad()
          + patch().deltaCoeffs()*
            (
                thermo.he(Tw, patchi)
              - thermo.he(pw, Tw, patch().faceCells())
            );
    }
    else if (isA<mixedEnergyCalculatedTemperatureFvPatchScalarField>(Tp))
    {
        mixedEnergyCalculatedTemperatureFvPatchScalarField& Tm =
            refCast<mixedEnergyCalculatedTemperatureFvPatchScalarField>(Tp);

        Tm.evaluate();

        valueFraction() = Tm.heValueFraction();
        refValue() = Tm.heRefValue();
        refGrad() = Tm.heRefGrad();
    }
    else
    {
        FatalErrorInFunction
            << "Temperature boundary condition not recognised."
            << "A " << typeName << " condition for energy must be used with a "
            << mixedFvPatchScalarField::typeName << " or "
            << mixedEnergyCalculatedTemperatureFvPatchScalarField::typeName
            << " condition for temperature."
            << exit(FatalError);
    }

    mixedFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        mixedEnergyFvPatchScalarField
    );
}

// ************************************************************************* //
