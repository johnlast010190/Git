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
    (c) 2017-2020 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "surfaceFilmModels/submodels/kinematic/force/contactAngleForces/distribution/distributionContactAngleForce.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(distributionContactAngleForce, 0);
addToRunTimeSelectionTable(force, distributionContactAngleForce, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

distributionContactAngleForce::distributionContactAngleForce
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    contactAngleForce(typeName, film, dict),
    rndGen_(),
    distribution_
    (
        distributionModel::New
        (
            coeffDict_.subDict("distribution"),
            rndGen_
        )
    ),
    theta_
    (
        volScalarField::New
        (
            IOobject::modelName("theta", typeName),
            filmModel_.regionMesh(),
            dimensionedScalar(dimless, 0)
        )
    ),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

distributionContactAngleForce::~distributionContactAngleForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> distributionContactAngleForce::theta() const
{
    if (curTimeIndex_ != filmModel_.time().timeIndex())
    {
        volScalarField::Internal& thetai = theta_;

        forAll(thetai, celli)
        {
            thetai[celli] = distribution_->sample();
        }

        forAll(theta_.boundaryField(), patchi)
        {
            if (!filmModel_.isCoupledPatch(patchi))
            {
                scalarField& thetap = theta_.boundaryFieldRef()[patchi];

                forAll(thetap, facei)
                {
                    thetap[facei] = distribution_->sample();
                }
            }
        }

        curTimeIndex_ = filmModel_.time().timeIndex();
    }

    return theta_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
