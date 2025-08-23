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

\*---------------------------------------------------------------------------*/

#include "relativeVelocityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(relativeVelocityModel, 0);
    defineRunTimeSelectionTable(relativeVelocityModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModel::relativeVelocityModel
(
    const dictionary& dict,
    const immiscibleIncompressibleTwoPhaseMixture& mixture
)
:
    mixture_(mixture),
    alphac_(mixture.alpha2()),
    alphad_(mixture.alpha1()),
    rhoc_(mixture.rho2()),
    rhod_(mixture.rho1()),
    Urel_
    (
        IOobject
        (
            "Urel",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector("Urel", dimVelocity, Zero)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::relativeVelocityModel> Foam::relativeVelocityModel::New
(
    const dictionary& dict,
    const immiscibleIncompressibleTwoPhaseMixture& mixture
)
{
    word modelType(dict.lookup(typeName));

    Info<< "Selecting relative velocity model " << modelType << endl;

    auto ctor = ctorTableLookup("time scale model", dictionaryConstructorTable_(), modelType);
    return
        autoPtr<relativeVelocityModel>
        (
            ctor
            (
                dict.subDict(modelType + "Coeffs"),
                mixture
            )
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModel::~relativeVelocityModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volSymmTensorField> Foam::relativeVelocityModel::tauDm() const
{
    const volScalarField limitedAlphad
    (
        "limitedAlphad",
        min(max(alphad_, scalar(0)), scalar(1))
    );

    volScalarField beta(limitedAlphad * (1.-limitedAlphad) * rhod_ * rhoc_ / mixture_.rho());

    return volSymmTensorField::New("tauDm", beta*sqr(Urel_));
}


// ************************************************************************* //
