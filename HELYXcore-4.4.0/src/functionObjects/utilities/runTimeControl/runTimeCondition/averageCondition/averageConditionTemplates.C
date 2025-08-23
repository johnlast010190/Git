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
    (c) 2015 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2020 Engys Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::runTimeControls::averageCondition::calc
(
    const word& fieldName,
    const scalar alpha,
    const scalar beta,
    bool& satisfied,
    bool& processed
)
{
    const word valueType =
        state_.objectResultType(functionObjectName_, fieldName);

    if (pTraits<Type>::typeName != valueType)
    {
        return;
    }

    Type currentValue =
        state_.getObjectResult<Type>(functionObjectName_, fieldName);

    const word meanName(fieldName + "Mean");

    Type meanValue = getConditionResult<Type>(meanName);
    meanValue = alpha*meanValue + beta*currentValue;

    scalar delta = mag(meanValue - currentValue);

    Log << "        " << meanName << ": " << meanValue
            << ", delta: " << delta << nl;

    if (Pstream::master() || !Pstream::parRun())
    {
        file() << tab << meanValue << tab << delta;
    }

    setConditionResult(meanName, meanValue);

    if (delta > tolerance_)
    {
        satisfied = false;
    }

    processed = true;
}


// ************************************************************************* //
