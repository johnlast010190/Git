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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2017-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "solidBodyMotionFunctions/rotatingStepMotion/rotatingStepMotion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/unitConversion/unitConversion.H"
#include "fields/volFields/volFields.H"
#include "fvMesh/fvMesh.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(rotatingStepMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        rotatingStepMotion,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        rotatingStepMotion,
        registry
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::rotatingStepMotion::rotatingStepMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime,
    const word& frameName
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime, frameName),
    origin_(SBMFCoeffs_.lookup("origin")),
    axis_(SBMFCoeffs_.lookup("axis")),
    theta_(SBMFCoeffs_.lookup<scalar>("theta")),
    theta0_(SBMFCoeffs_.lookupOrDefault<scalar>("theta0", 0.0)),
    stepPeriod_(SBMFCoeffs_.lookup<scalar>("period")),
    nStepsDone_(-1)
{}


Foam::solidBodyMotionFunctions::rotatingStepMotion::rotatingStepMotion
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    solidBodyMotionFunction(obr, SBMFCoeffs, frameName),
    origin_(CofR0()),
    axis_(frame().axis0()),
    theta_(SBMFCoeffs_.lookup<scalar>("theta")),
    theta0_(SBMFCoeffs_.lookupOrDefault<scalar>("theta0", 0.0)),
    stepPeriod_(SBMFCoeffs_.lookup<scalar>("period")),
    nStepsDone_(-1)
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::rotatingStepMotion::~rotatingStepMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::rotatingStepMotion::transformation
(
    const scalar t1,
    const scalar t2
) const
{
    // Index counting should be the same as time but can be generalised and used
    // for incremental motion
    label nSteps(t2/stepPeriod_);

    if (hasFrame() && isIncrementalMotion())
    {
        theta0_ = (nStepsDone_ == -1) ? theta0_ : 0.0;
        if (nStepsDone_ < nSteps)
        {
            nSteps -= nStepsDone_;
            nStepsDone_ += nSteps;
        }

        origin_ = frame().coorSys().origin();
        axis_ = normalised(frame().axis());
    }

    // Stepwise rotation around axis
    const scalar angle = degToRad(theta0_ + nSteps*theta_);

    const quaternion R(axis_, angle);
    const septernion TR = (septernion(-origin_)*R*septernion(origin_));

    return TR;
}


Foam::vectorTuple
Foam::solidBodyMotionFunctions::rotatingStepMotion::velocity() const
{
    // Dummy steady-state motion, so return zero velocity!
    return vectorTuple(Zero, Zero);
}


bool Foam::solidBodyMotionFunctions::rotatingStepMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    return true;
}

// ************************************************************************* //
