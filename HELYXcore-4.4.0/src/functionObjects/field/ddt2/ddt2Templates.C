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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "db/Time/Time.H"
#include "dimensionedTypes/dimensionedType/dimensionedType.H"
#include "finiteVolume/fvc/fvcDdt.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FieldType>
int Foam::functionObjects::ddt2::apply(const word& inputName, int& state)
{
    // State: return 0 (not-processed), -1 (skip), +1 ok

    // Already done, or not available
    if (state || !foundObject<FieldType>(inputName))
    {
        return state;
    }

    const FieldType& input = lookupObject<FieldType>(inputName);

    word outputName(resultName_);
    outputName.replace("@@", inputName);

    results_.set(outputName);

    if (!foundObject<volScalarField>(outputName))
    {
        tmp<volScalarField> tddt2
        (
            new volScalarField
            (
                IOobject
                (
                    outputName,
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    (
                        mag_
                      ? Foam::mag(input.dimensions()/dimTime)
                      : Foam::magSqr(input.dimensions()/dimTime)
                    ),
                    0
                )
            )
        );

        store(outputName, tddt2);
    }

    volScalarField& output =
        const_cast<volScalarField&>(lookupObject<volScalarField>(outputName));

    if (mag_)
    {
        output = Foam::mag(fvc::ddt(input));
    }
    else
    {
        output = Foam::magSqr(fvc::ddt(input));
    }

    // Could add additional statistics here
    Log << type() << ' ' << name()
        << " field " << outputName
        << " average: " << gAverage(output) << endl;

    state = +1;
    return state;
}


// ************************************************************************* //
