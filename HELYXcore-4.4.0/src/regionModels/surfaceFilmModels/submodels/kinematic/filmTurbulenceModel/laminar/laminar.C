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
    (c) 2011-2021 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "surfaceFilmModels/submodels/kinematic/filmTurbulenceModel/laminar/laminar.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "db/Time/Time.H"
#include "fields/volFields/volFields.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "surfaceFilmModels/kinematicSingleLayer/kinematicSingleLayer.H"
#include "fields/fvPatchFields/basic/extrapolatedCalculated/extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(laminar, 0);
addToRunTimeSelectionTable(turbulenceModel, laminar, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

laminar::laminar
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    turbulenceModel(type(), film, dict),
    Cf_(coeffDict_.lookup<scalar>("Cf"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

laminar::~laminar()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volVectorField::Internal> laminar::Us() const
{
    tmp<volVectorField::Internal> tUs
    (
        volVectorField::Internal::New
        (
            IOobject::modelName("Us", typeName),
            1.5*filmModel_.U()
        )
    );

    return tUs;
}


void laminar::correct()
{}


tmp<fvVectorMatrix> laminar::Su(volVectorField& U) const
{
    // Local reference to film model
    const kinematicSingleLayer& film =
        static_cast<const kinematicSingleLayer&>(filmModel_);

    // Local references to film fields
    const volScalarField::Internal& mu = film.mu();
    const volScalarField::Internal& rho = film.rho();
    const volScalarField::Internal& delta = film.delta();
    const volVectorField::Internal& Up = film.UPrimary();
    const volScalarField::Internal& rhop = film.rhoPrimary();
    const volScalarField::Internal& VbyA = film.VbyA();

    // Employ simple coeff-based model
    volScalarField::Internal Cs("Cs", Cf_*rhop*mag(Up - U)/VbyA);
    volScalarField::Internal Cw
    (
        "Cw",
        mu/((1.0/3.0)*VbyA*delta + mu*film.time().deltaT()/rho)
    );

    return
    (
       - fvm::Sp(Cs, U) + Cs*Up // surface contribution
       - fvm::Sp(Cw, U) + Cw*film.Uw() // wall contribution
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
