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
    (c) 2024 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/mappedField/mappedPatchFieldBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void mappedPatchFieldBase<vector>::scaleWithAverage(Field<vector>& fieldValue) const
{
    if (setAverage_ && mag(flowRate_)<SMALL)
    {
        vector averagePsi =
            gSum(patchField_.patch().magSf()*fieldValue)
           /gSum(patchField_.patch().magSf());

        if (mag(averagePsi)/mag(average_) > 0.5)
        {
            fieldValue *= mag(average_)/mag(averagePsi);
        }
        else
        {
            fieldValue += (average_ - averagePsi);
        }
    }

    else if (setAverage_ && mag(flowRate_)>SMALL)
    {
        //- target U from target flowRate
        scalar Umean = getUfromFlowRate();

        tmp<vectorField> nf(patchField_.patch().nf());
        vectorField Un((fieldValue&patchField_.patch().nf())*patchField_.patch().nf());
        vectorField Ut(fieldValue -Un);

        //current flowRate
        scalar currentFlowRate;

        // current Umean from current flowRate
        scalar averagePsi;

        if (isMassFlow_)
        {
            const fvPatchField <scalar> *rhoP =
                &patchField_.patch().lookupPatchFieldInDb<volScalarField, scalar>
                (
                    patchField_.db(),
                    "rho"
                );

            scalar rhoMean =
                gSum((*rhoP)*(*rhoP).patch().magSf())/gSum((*rhoP).patch().magSf());

            currentFlowRate =
                gSum((patchField_.patch().Sf() & fieldValue)*rhoMean);

            averagePsi = currentFlowRate/gSum(patchField_.patch().magSf())/rhoMean;
        }
        else
        {
            currentFlowRate =
                gSum(patchField_.patch().Sf() & fieldValue);

            averagePsi = currentFlowRate/gSum(patchField_.patch().magSf());
        }

        if (mag(averagePsi)/mag(Umean) > 0.5)
        {
            Un *= mag(Umean)/mag(averagePsi);
        }
        else
        {
            Un += ((Umean*nf()) - (averagePsi*nf()));
        }
        fieldValue = Un +Ut;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
