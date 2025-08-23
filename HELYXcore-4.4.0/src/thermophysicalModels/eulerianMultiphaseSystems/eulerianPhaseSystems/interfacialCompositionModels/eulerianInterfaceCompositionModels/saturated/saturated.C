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

#include "saturated.H"
#include "eulerianPhaseSystems/eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianInterfaceCompositionModels
{
    defineTypeNameAndDebug(saturated, 0);
    addToRunTimeSelectionTable
    (
        eulerianInterfaceCompositionModel,
        saturated,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::eulerianInterfaceCompositionModels::saturated::wRatioByP() const
{
    return
        composition().W(saturatedIndex_)/composition().W()/thermo().p();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianInterfaceCompositionModels::saturated::saturated
(
    const dictionary& dict,
    const eulerianPhasePair& pair
)
:
    eulerianInterfaceCompositionModel(dict, pair),
    saturatedName_(species()[0]),
    saturatedIndex_(composition().species()[saturatedName_]),
    saturationModel_
    (
        eulerianSaturationModel::New(dict.subDict("saturationPressure"))
    )
{
    if (species().size() != 1)
    {
        FatalErrorInFunction
            << "saturated model is suitable for one species only."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianInterfaceCompositionModels::saturated::~saturated()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::eulerianInterfaceCompositionModels::saturated::update
(
    const volScalarField& Tf
)
{}


Foam::tmp<Foam::volScalarField>
Foam::eulerianInterfaceCompositionModels::saturated::Yf
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    if (saturatedName_ == speciesName)
    {
        return wRatioByP()*saturationModel_->pSat(Tf);
    }
    else
    {
        const label speciesIndex(composition().species()[speciesName]);

        return
            composition().Y()[speciesIndex]
           *(scalar(1) - wRatioByP()*saturationModel_->pSat(Tf))
           /max(scalar(1) - composition().Y()[saturatedIndex_], SMALL);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::eulerianInterfaceCompositionModels::saturated::YfPrime
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    if (saturatedName_ == speciesName)
    {
        return wRatioByP()*saturationModel_->pSatPrime(Tf);
    }
    else
    {
        const label speciesIndex = composition().species()[speciesName];

        return
          - composition().Y()[speciesIndex]
           *wRatioByP()*saturationModel_->pSatPrime(Tf)
           /max(scalar(1) - composition().Y()[saturatedIndex_], SMALL);
    }
}


// ************************************************************************* //
