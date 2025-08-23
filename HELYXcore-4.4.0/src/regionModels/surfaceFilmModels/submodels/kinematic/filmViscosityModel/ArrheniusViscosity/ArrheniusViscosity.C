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
    (c) 2015-2020 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "surfaceFilmModels/submodels/kinematic/filmViscosityModel/ArrheniusViscosity/ArrheniusViscosity.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ArrheniusViscosity, 0);
addToRunTimeSelectionTable(viscosityModel, ArrheniusViscosity, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ArrheniusViscosity::ArrheniusViscosity
(
    surfaceFilmRegionModel& film,
    const dictionary& dict,
    volScalarField& mu
)
:
    viscosityModel(typeName, film, dict, mu),
    viscosity_(viscosityModel::New(film, coeffDict_, mu)),
    k1_("k1", dimTemperature, coeffDict_),
    k2_("k2", dimTemperature, coeffDict_),
    Tref_("Tref", dimTemperature, coeffDict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ArrheniusViscosity::~ArrheniusViscosity()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void ArrheniusViscosity::correct
(
    const volScalarField& p,
    const volScalarField& T
)
{
    viscosity_->correct(p, T);
    mu_ *= exp(k1_*((1/(T + k2_)) - 1/(Tref_ + k2_)));
    mu_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
