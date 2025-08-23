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
    (c) 2017-2021 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "surfaceFilmModels/submodels/kinematic/force/contactAngleForces/temperatureDependent/temperatureDependentContactAngleForce.H"
#include "surfaceFilmModels/thermoSingleLayer/thermoSingleLayer.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(temperatureDependentContactAngleForce, 0);
addToRunTimeSelectionTable
(
    force,
    temperatureDependentContactAngleForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

temperatureDependentContactAngleForce::temperatureDependentContactAngleForce
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    contactAngleForce(typeName, film, dict),
    thetaPtr_(Function1<scalar>::New("theta", coeffDict_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

temperatureDependentContactAngleForce::~temperatureDependentContactAngleForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> temperatureDependentContactAngleForce::theta() const
{
    tmp<volScalarField> ttheta
    (
        volScalarField::New
        (
            IOobject::modelName("theta", typeName),
            filmModel_.regionMesh(),
            dimensionedScalar(dimless, 0)
        )
    );

    volScalarField& theta = ttheta.ref();

    const thermoSingleLayer& film = filmType<thermoSingleLayer>();

    const volScalarField& T = film.thermo().T();

    theta.ref().field() = thetaPtr_->value(T());

    forAll(theta.boundaryField(), patchi)
    {
        if (!filmModel_.isCoupledPatch(patchi))
        {
            theta.boundaryFieldRef()[patchi] =
                thetaPtr_->value(T.boundaryField()[patchi]);
        }
    }

    return ttheta;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
