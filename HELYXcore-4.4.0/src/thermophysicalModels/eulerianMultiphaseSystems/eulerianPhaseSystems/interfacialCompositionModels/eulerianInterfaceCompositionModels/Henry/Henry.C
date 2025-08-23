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

#include "Henry.H"
#include "eulerianPhaseSystems/eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianInterfaceCompositionModels
{
    defineTypeNameAndDebug(Henry, 0);
    addToRunTimeSelectionTable
    (
        eulerianInterfaceCompositionModel,
        Henry,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianInterfaceCompositionModels::Henry::Henry
(
    const dictionary& dict,
    const eulerianPhasePair& pair
)
:
    eulerianInterfaceCompositionModel(dict, pair),
    k_(dict.lookup("k")),
    YSolvent_
    (
        IOobject
        (
            IOobject::groupName("YSolvent", pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh()
        ),
        pair.phase1().mesh(),
        dimensionedScalar(dimless, 1)
    )
{
    if (k_.size() != species().size())
    {
        FatalErrorInFunction
            << "Differing number of species and solubilities"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianInterfaceCompositionModels::Henry::~Henry()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::eulerianInterfaceCompositionModels::Henry::update
(
    const volScalarField& Tf
)
{
    YSolvent_ = scalar(1);

    forAllConstIter(hashedWordList, species(), iter)
    {
        YSolvent_ -= Yf(*iter, Tf);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::eulerianInterfaceCompositionModels::Henry::Yf
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    if (species().found(speciesName))
    {
        const label index = species()[speciesName];

        return
            k_[index]
           *otherComposition().Y(speciesName)
           *otherThermo().rho()
           /thermo().rho();
    }
    else
    {
        return YSolvent_*composition().Y(speciesName);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::eulerianInterfaceCompositionModels::Henry::YfPrime
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    return volScalarField::New
    (
        IOobject::groupName("YfPrime", pair().name()),
        pair().phase1().mesh(),
        dimensionedScalar(dimless/dimTemperature, 0)
    );
}


// ************************************************************************* //
