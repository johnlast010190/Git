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
    (c) 2010-2012 Engys Ltd.
    (c) 2011-2020 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "turbulentFluidThermoModels/derivedFvPatchFields/wallFunctions/DtWallFunctions/DtWallFunction/DtWallFunctionFvPatchScalarField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DtWallFunctionFvPatchScalarField::
DtWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


DtWallFunctionFvPatchScalarField::
DtWallFunctionFvPatchScalarField
(
    const DtWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


DtWallFunctionFvPatchScalarField::
DtWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{}


DtWallFunctionFvPatchScalarField::
DtWallFunctionFvPatchScalarField
(
    const DtWallFunctionFvPatchScalarField& awfpsf
)
:
    fixedValueFvPatchScalarField(awfpsf)
{}


DtWallFunctionFvPatchScalarField::
DtWallFunctionFvPatchScalarField
(
    const DtWallFunctionFvPatchScalarField& awfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(awfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void DtWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Retrieve turbulence properties from model
    const compressible::turbulenceModel& turbModel =
        db().lookupObject<compressible::turbulenceModel>
        (
            IOobject::groupName
            (
                compressible::turbulenceModel::propertiesName,
                internalField().group()
            )
        );

    tmp<scalarField> tmutw = turbModel.mut(patch().index());
    const scalarField& mutw = tmutw();

    string pofix(string(this->internalField().name())(2, 100));
    const word subDictName(pofix + "concentrationTransport");
    bool isMat = (basicThermo::dictName == basicThermo::matDictName);

    //read turbulent Schmidt number
    dimensionedScalar Sct
    (
        isMat ? "Sct" : "Sct" + pofix,
        isMat
      ? turbModel.transport().properties().optionalSubDict(subDictName).lookup("Sct")
      : turbModel.transport().properties().lookup("Sct" + pofix)
    );

    forceAssign(mutw/Sct.value());

    fixedValueFvPatchScalarField::updateCoeffs();

}


void DtWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, DtWallFunctionFvPatchScalarField);

// Backward compatibility
addSpecialNamedToPatchFieldRunTimeSelection
(
    fvPatchScalarField,
    DtWallFunctionFvPatchScalarField,
    compressible::DtWallFunction,
    compressible__DtWallFunction
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
