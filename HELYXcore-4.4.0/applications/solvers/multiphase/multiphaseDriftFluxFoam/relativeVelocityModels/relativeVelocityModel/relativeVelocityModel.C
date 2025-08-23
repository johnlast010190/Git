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
    (c) 2014-2019 OpenFOAM Foundation
    (c) 2020 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "relativeVelocityModel.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "fields/fvPatchFields/derived/slip/slipFvPatchFields.H"
#include "fields/fvPatchFields/derived/partialSlip/partialSlipFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(relativeVelocityModel, 0);
    defineRunTimeSelectionTable(relativeVelocityModel, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions   * * * * * * * * * * * //

Foam::wordList Foam::relativeVelocityModel::UdmPatchFieldTypes() const
{
    const volVectorField& U = alpha_.mesh().lookupObject<volVectorField>("U");

    wordList UdmTypes
    (
        U.boundaryField().size(),
        calculatedFvPatchScalarField::typeName
    );

    forAll(U.boundaryField(), i)
    {
        if
        (
            isA<fixedValueFvPatchVectorField>(U.boundaryField()[i])
         || isA<slipFvPatchVectorField>(U.boundaryField()[i])
         || isA<partialSlipFvPatchVectorField>(U.boundaryField()[i])
        )
        {
            UdmTypes[i] = fixedValueFvPatchVectorField::typeName;
        }
    }

    return UdmTypes;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModel::relativeVelocityModel
(
    const word& phaseName,
    const dictionary& dict,
    const volScalarField& alpha,
    const scalar& alphaMax
)
:
    alpha_(alpha),
    alphaMax_(alphaMax),
    Udm_
    (
        IOobject
        (
            IOobject::groupName("Udm", phaseName),
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha_.mesh(),
        dimensionedVector("Udm", dimVelocity, Zero),
        UdmPatchFieldTypes()
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::relativeVelocityModel> Foam::relativeVelocityModel::New
(
    const word& phaseName,
    const dictionary& dict,
    const volScalarField& alpha,
    const scalar& alphaMax
)
{
    const word modelType(dict.lookup(typeName));

    Info<< "Selecting relative velocity model " << modelType << endl;

    const auto ctor =
        ctorTableLookup
        (
            "time scale model type",
            dictionaryConstructorTable_(),
            modelType
        );
    return
        autoPtr<relativeVelocityModel>
        (
            ctor
            (
                phaseName,
                dict.optionalSubDict(modelType + "Coeffs"),
                alpha,
                alphaMax
            )
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModel::~relativeVelocityModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
