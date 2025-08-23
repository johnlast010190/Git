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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/acousticWaveTransmissive/acousticWaveTransmissiveFvPatchField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "finiteVolume/ddtSchemes/EulerDdtScheme/EulerDdtScheme.H"
#include "finiteVolume/ddtSchemes/CrankNicolsonDdtScheme/CrankNicolsonDdtScheme.H"
#include "finiteVolume/ddtSchemes/backwardDdtScheme/backwardDdtScheme.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::acousticWaveTransmissiveFvPatchField<Type>::acousticWaveTransmissiveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(p, iF),
    advectiveU_(0)
{}


template<class Type>
Foam::acousticWaveTransmissiveFvPatchField<Type>::acousticWaveTransmissiveFvPatchField
(
    const acousticWaveTransmissiveFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    advectiveFvPatchField<Type>(ptf, p, iF, mapper),
    advectiveU_(ptf.advectiveU_)
{}


template<class Type>
Foam::acousticWaveTransmissiveFvPatchField<Type>::acousticWaveTransmissiveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    advectiveFvPatchField<Type>(p, iF, dict),
    advectiveU_(dict.lookup<scalar>("advectiveSpeed"))
{}


template<class Type>
Foam::acousticWaveTransmissiveFvPatchField<Type>::acousticWaveTransmissiveFvPatchField
(
    const acousticWaveTransmissiveFvPatchField& ptpsf
)
:
    advectiveFvPatchField<Type>(ptpsf),
    advectiveU_(ptpsf.advectiveU_)
{}


template<class Type>
Foam::acousticWaveTransmissiveFvPatchField<Type>::acousticWaveTransmissiveFvPatchField
(
    const acousticWaveTransmissiveFvPatchField& ptpsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(ptpsf, iF),
    advectiveU_(ptpsf.advectiveU_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::scalarField>
Foam::acousticWaveTransmissiveFvPatchField<Type>::advectionSpeed() const
{
    return tmp<scalarField>::New(this->size(), advectiveU_);
}

// template<class Type>
// void Foam::acousticWaveTransmissiveFvPatchField<Type>::updateCoeffs()
// {
//     Info<< "DIO " << endl;
//     if
//     (
//         this->db().objectRegistry::template foundObject<surfaceScalarField>
//         (advectiveFvPatchField<Type>::phiName_)
//      && this->db().objectRegistry::template foundObject<volScalarField>
//         (advectiveFvPatchField<Type>::rhoName_)
//      && this->db().objectRegistry::template foundObject<volScalarField>
//         (psiName_)
//     )
//     {
//         if (totalPressure_) totalPressure();

//         advectiveFvPatchField<Type>::updateCoeffs();
//     }
//     else
//     {
//         this->valueFraction() = 1.0;
//         mixedFvPatchField<Type>::updateCoeffs();
//     }
// }

template<class Type>
void Foam::acousticWaveTransmissiveFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeEntry("advectiveSpeed", advectiveU_);
    // fvPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
