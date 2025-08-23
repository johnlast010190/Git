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
    (c) 2015-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "interfaceCompositionModel.H"
#include "phaseModel/phaseModel/phaseModel.H"
#include "phasePair/phasePair/phasePair.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interfaceCompositionModel, 0);
    defineRunTimeSelectionTable(interfaceCompositionModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCompositionModel::interfaceCompositionModel
(
    const dictionary& dict,
    const phasePair& pair
)
:
    pair_(pair),
    species_(dict.lookup("species")),
    Le_("Le", dimless, dict),
    thermo_
    (
        pair.phase1().mesh().lookupObject<rhoMulticomponentThermo>
        (
            IOobject::groupName(basicThermo::dictName, pair.phase1().name())
        )
    ),
    otherThermo_
    (
        pair.phase2().mesh().lookupObject<rhoThermo>
        (
            IOobject::groupName(basicThermo::dictName, pair.phase2().name())
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceCompositionModel::~interfaceCompositionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::interfaceCompositionModel::dY
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    const label speciei = composition().species()[speciesName];

    return Yf(speciesName, Tf) - composition().Y()[speciei];
}


Foam::tmp<Foam::volScalarField> Foam::interfaceCompositionModel::D
(
    const word& speciesName
) const
{
    const label speciei = composition().species()[speciesName];
    const volScalarField& p = thermo_.p();
    const volScalarField& T = thermo_.T();

    return volScalarField::New
    (
        IOobject::groupName("D", pair_.name()),
        composition().kappa(speciei, p, T)
       /composition().Cp(speciei, p, T)
       /composition().rho(speciei, p, T)
       /Le_
    );
}


Foam::tmp<Foam::volScalarField> Foam::interfaceCompositionModel::L
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    const label speciei = composition().species()[speciesName];
    const volScalarField& p = thermo_.p();
    volScalarField Ha(composition().ha(speciei, p, Tf));

    const volScalarField& otherP(otherThermo_.p());
    tmp<volScalarField> otherHa(nullptr);
    if (otherHasComposition())
    {
        const label otherSpeciei = otherComposition().species()[speciesName];
        otherHa = otherComposition().ha(otherSpeciei, otherP, Tf);
    }
    else
    {
        otherHa = otherThermo_.ha(otherP, Tf);
    }

    return
        volScalarField::New
        (
            IOobject::groupName("L", pair_.name()),
            Ha - otherHa
        );
}


void Foam::interfaceCompositionModel::addMDotL
(
    const volScalarField& K,
    const volScalarField& Tf,
    volScalarField& mDotL,
    volScalarField& mDotLPrime
) const
{
    forAllConstIter(hashedWordList, species_, iter)
    {
        const volScalarField rhoKDL(thermo_.rho()*K*D(*iter)*L(*iter, Tf));

        mDotL += rhoKDL*dY(*iter, Tf);
        mDotLPrime += rhoKDL*YfPrime(*iter, Tf);
    }
}


// ************************************************************************* //
