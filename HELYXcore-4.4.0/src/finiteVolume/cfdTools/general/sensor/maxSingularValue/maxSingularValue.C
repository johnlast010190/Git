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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2015 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/sensor/maxSingularValue/maxSingularValue.H"
#include "cfdTools/general/include/fvCFD.H"
#include "global/constants/mathematical/mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::sensorTypes::maxSingularValue<Type>::maxSingularValue
(
    const fvMesh& mesh,
    const dictionary& sensorDict
)
:
    sensor<Type>(mesh, sensorDict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::sensorTypes::maxSingularValue<Type>::~maxSingularValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Foam::scalar>>
Foam::sensorTypes::maxSingularValue<Type>::valueField() const
{
    const volVectorField& U =
        this->mesh_.objectRegistry::template lookupObject<volVectorField>
        (
            this->fieldName()
        );

    tensorField operatorT = fvc::grad(U)().primitiveField();
    operatorT = operatorT.T() & operatorT;

    vectorField eigenValuesField ( eigenValues(operatorT) );

    return tmp<scalarField>
    (
        new scalarField(sqrt(eigenValuesField.component(2)))
    );
}


// ************************************************************************* //
