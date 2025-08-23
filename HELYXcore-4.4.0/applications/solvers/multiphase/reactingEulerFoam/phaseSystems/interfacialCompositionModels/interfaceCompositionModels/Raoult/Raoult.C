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

#include "Raoult.H"
#include "phasePair/phasePair/phasePair.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interfaceCompositionModels
{
    defineTypeNameAndDebug(Raoult, 0);
    addToRunTimeSelectionTable(interfaceCompositionModel, Raoult, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCompositionModels::Raoult::Raoult
(
    const dictionary& dict,
    const phasePair& pair
)
:
    interfaceCompositionModel(dict, pair),
    YNonVapour_
    (
        IOobject
        (
            IOobject::groupName("YNonVapour", pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh()
        ),
        pair.phase1().mesh(),
        dimensionedScalar(dimless, 1)
    ),
    YNonVapourPrime_
    (
        IOobject
        (
            IOobject::groupName("YNonVapourPrime", pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh()
        ),
        pair.phase1().mesh(),
        dimensionedScalar(dimless/dimTemperature, 0)
    )
{
    forAllConstIter(hashedWordList, species(), iter)
    {
        speciesModels_.insert
        (
            *iter,
            autoPtr<interfaceCompositionModel>
            (
                interfaceCompositionModel::New(dict.subDict(*iter), pair)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceCompositionModels::Raoult::~Raoult()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::interfaceCompositionModels::Raoult::update(const volScalarField& Tf)
{
    YNonVapour_ = scalar(1);

    forAllIter
    (
        HashTable<autoPtr<interfaceCompositionModel>>,
        speciesModels_,
        iter
    )
    {
        iter()->update(Tf);
        YNonVapour_ -=
            otherComposition().Y(iter.key())*iter()->Yf(iter.key(), Tf);
        YNonVapourPrime_ -=
            otherComposition().Y(iter.key())*iter()->YfPrime(iter.key(), Tf);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceCompositionModels::Raoult::Yf
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    if (species().found(speciesName))
    {
        return
            otherComposition().Y(speciesName)
           *speciesModels_[speciesName]->Yf(speciesName, Tf);
    }
    else
    {
        return composition().Y(speciesName)*YNonVapour_;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceCompositionModels::Raoult::YfPrime
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    if (species().found(speciesName))
    {
        return
            otherComposition().Y(speciesName)
           *speciesModels_[speciesName]->YfPrime(speciesName, Tf);
    }
    else
    {
        return
            otherComposition().Y(speciesName)
           *YNonVapourPrime_;
    }
}


// ************************************************************************* //
