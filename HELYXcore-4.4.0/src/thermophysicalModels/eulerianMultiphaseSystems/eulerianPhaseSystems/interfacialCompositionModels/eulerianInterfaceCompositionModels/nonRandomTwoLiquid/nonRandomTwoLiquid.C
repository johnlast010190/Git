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
    (c) 2015-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "nonRandomTwoLiquid.H"
#include "eulerianPhaseSystems/eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianInterfaceCompositionModels
{
    defineTypeNameAndDebug(nonRandomTwoLiquid, 0);
    addToRunTimeSelectionTable
    (
        eulerianInterfaceCompositionModel,
        nonRandomTwoLiquid,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianInterfaceCompositionModels::nonRandomTwoLiquid::nonRandomTwoLiquid
(
    const dictionary& dict,
    const eulerianPhasePair& pair
)
:
    eulerianInterfaceCompositionModel(dict, pair),
    gamma1_
    (
        IOobject
        (
            IOobject::groupName("gamma1", pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh()
        ),
        pair.phase1().mesh(),
        dimensionedScalar("one", dimless, 1)
    ),
    gamma2_
    (
        IOobject
        (
            IOobject::groupName("gamma2", pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh()
        ),
        pair.phase1().mesh(),
        dimensionedScalar("one", dimless, 1)
    ),
    beta12_("", dimless/dimTemperature, 0),
    beta21_("", dimless/dimTemperature, 0)
{
    if (species().size() != 2)
    {
        FatalErrorInFunction
            << "nonRandomTwoLiquid model is suitable for two species only."
            << exit(FatalError);
    }

    species1Name_ = species()[0];
    species2Name_ = species()[1];

    species1Index_ = composition().species()[species1Name_];
    species2Index_ = composition().species()[species2Name_];

    alpha12_ = dimensionedScalar
    (
        "alpha12",
        dimless,
        dict.subDict(species1Name_).lookup("alpha")
    );
    alpha21_ = dimensionedScalar
    (
        "alpha21",
        dimless,
        dict.subDict(species2Name_).lookup("alpha")
    );

    beta12_ = dimensionedScalar
    (
        "beta12",
        dimless/dimTemperature,
        dict.subDict(species1Name_).lookup("beta")
    );
    beta21_ = dimensionedScalar
    (
        "beta21",
        dimless/dimTemperature,
        dict.subDict(species2Name_).lookup("beta")
    );

    saturationModel12_.reset
    (
        eulerianSaturationModel::New
        (
            dict.subDict(species1Name_).subDict("interaction")
        ).ptr()
    );
    saturationModel21_.reset
    (
        eulerianSaturationModel::New
        (
            dict.subDict(species2Name_).subDict("interaction")
        ).ptr()
    );

    speciesModel1_.reset
    (
        eulerianInterfaceCompositionModel::New
        (
            eulerianInterfaceCompositionModel::typeName,
            dict.subDict(species1Name_),
            pair
        ).ptr()
    );
    speciesModel2_.reset
    (
        eulerianInterfaceCompositionModel::New
        (
            eulerianInterfaceCompositionModel::typeName,
            dict.subDict(species2Name_),
            pair
        ).ptr()
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianInterfaceCompositionModels::nonRandomTwoLiquid::
~nonRandomTwoLiquid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::eulerianInterfaceCompositionModels::nonRandomTwoLiquid::update
(
    const volScalarField& Tf
)
{
    const volScalarField W(thermo().composition().W());

    const volScalarField X1
    (
        composition().Y(species1Index_)*W/composition().W(species1Index_)
    );

    const volScalarField X2
    (
        composition().Y(species2Index_)*W/composition().W(species2Index_)
    );

    const volScalarField alpha12(alpha12_ + Tf*beta12_);
    const volScalarField alpha21(alpha21_ + Tf*beta21_);

    const volScalarField tau12(saturationModel12_->lnPSat(Tf));
    const volScalarField tau21(saturationModel21_->lnPSat(Tf));

    const volScalarField G12(exp(-alpha12*tau12));
    const volScalarField G21(exp(-alpha21*tau21));

    gamma1_ =
        exp
        (
            sqr(X2)
           *(
                tau21*sqr(G21)/max(sqr(X1 + X2*G21), SMALL)
              + tau12*G12/max(sqr(X2 + X1*G12), SMALL)
            )
        );
    gamma2_ =
        exp
        (
            sqr(X1)
           *(
                tau12*sqr(G12)/max(sqr(X2 + X1*G12), SMALL)
              + tau21*G21/max(sqr(X1 + X2*G21), SMALL)
            )
        );
}


Foam::tmp<Foam::volScalarField>
Foam::eulerianInterfaceCompositionModels::nonRandomTwoLiquid::Yf
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    if (speciesName == species1Name_)
    {
        return
            otherComposition().Y(speciesName)
           *speciesModel1_->Yf(speciesName, Tf)
           *gamma1_;
    }
    else if (speciesName == species2Name_)
    {
        return
            otherComposition().Y(speciesName)
           *speciesModel2_->Yf(speciesName, Tf)
           *gamma2_;
    }
    else
    {
        return
            composition().Y(speciesName)
           *(scalar(1) - Yf(species1Name_, Tf) - Yf(species2Name_, Tf));
    }
}


Foam::tmp<Foam::volScalarField>
Foam::eulerianInterfaceCompositionModels::nonRandomTwoLiquid::YfPrime
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    if (speciesName == species1Name_)
    {
        return
            otherComposition().Y(speciesName)
           *speciesModel1_->YfPrime(speciesName, Tf)
           *gamma1_;
    }
    else if (speciesName == species2Name_)
    {
        return
            otherComposition().Y(speciesName)
           *speciesModel2_->YfPrime(speciesName, Tf)
           *gamma2_;
    }
    else
    {
        return
          - composition().Y(speciesName)
           *(YfPrime(species1Name_, Tf) + YfPrime(species2Name_, Tf));
    }
}


// ************************************************************************* //
