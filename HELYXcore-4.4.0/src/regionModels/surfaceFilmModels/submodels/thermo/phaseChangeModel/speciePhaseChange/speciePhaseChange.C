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
    (c) 2021 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "surfaceFilmModels/submodels/thermo/phaseChangeModel/speciePhaseChange/speciePhaseChange.H"
#include "surfaceFilmModels/thermoSingleLayer/thermoSingleLayer.H"
#include "fluidThermo/fluidThermo.H"
#include "mixtures/basicSpecieMixture/basicSpecieMixture.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(speciePhaseChange, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

speciePhaseChange::speciePhaseChange
(
    const word& modelType,
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    phaseChangeModel(modelType, film, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

speciePhaseChange::~speciePhaseChange()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label speciePhaseChange::vapId() const
{
    const thermoSingleLayer& film = filmType<thermoSingleLayer>();

    // Set local liquidThermo properties
    const basicSpecieMixture& primarySpecieThermo =
        refCast<const basicSpecieMixture>(film.primaryThermo());
    const word evaporationSpecies =
        film.thermo().properties().lookup<word>("evaporationSpecies");
    return primarySpecieThermo.species()[evaporationSpecies];
}


Foam::scalar speciePhaseChange::Wvap() const
{
    const thermoSingleLayer& film = filmType<thermoSingleLayer>();

    const basicSpecieMixture& primarySpecieThermo =
        refCast<const basicSpecieMixture>(film.primaryThermo());

    return primarySpecieThermo.W(vapId());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //