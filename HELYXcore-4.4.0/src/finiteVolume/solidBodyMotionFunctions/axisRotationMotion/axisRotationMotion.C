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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2019-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "solidBodyMotionFunctions/axisRotationMotion/axisRotationMotion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/unitConversion/unitConversion.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(axisRotationMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        axisRotationMotion,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        axisRotationMotion,
        registry
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::axisRotationMotion::axisRotationMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime,
    const word& frameName
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime, frameName),
    origin_(SBMFCoeffs_.lookup("origin"))
{
    axisRotationMotion::read(SBMFCoeffs);
}


Foam::solidBodyMotionFunctions::axisRotationMotion::axisRotationMotion
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    solidBodyMotionFunction(obr, SBMFCoeffs, frameName),
    origin_(CofR0())
{
    axisRotationMotion::read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::axisRotationMotion::~axisRotationMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::axisRotationMotion::transformation
(
    const scalar t1,
    const scalar t2
) const
{
    const scalar dt = t2 - t1;

    septernion TR = septernion::I;

    if (mag(omega_) > VSMALL)
    {
        vector axis = normalised(omega_);

        if (hasFrame())
        {
            if (isIncrementalMotion())
            {
                axis = frame().cartesianSys().transform(axis);
                origin_ = frame().coorSys().origin();
            }
            else
            {
                axis = frame().cartesianSys0().transform(axis);
            }
        }

        scalar angle = mag(omega_*dt);

        if (angle > VSMALL)
        {
            const quaternion R(axis, angle);
            TR = septernion(-origin_)*R*septernion(origin_);
        }
    }

    DebugInFunction
        << "Transformation from time " << t1
        << " to " << t2 << " is " << TR << endl;

    return TR;
}


Foam::vectorTuple
Foam::solidBodyMotionFunctions::axisRotationMotion::velocity() const
{
    if (mag(omega_) > VSMALL)
    {
        vector axis(normalised(omega_));

        if (hasFrame())
        {
            axis = frame().cartesianSys().transform(normalised(omega_));
        }

        return vectorTuple(Zero, axis*mag(omega_));
    }
    else
    {
        return vectorTuple(Zero, Zero);
    }

}


bool Foam::solidBodyMotionFunctions::axisRotationMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    // Read rotational velocity in deg/s (should it be rad/s?)
    const vector rotVelocity(SBMFCoeffs_.lookup<vector>("radialVelocity"));
    omega_ =
        vector
        (
            degToRad(rotVelocity.x()),
            degToRad(rotVelocity.y()),
            degToRad(rotVelocity.z())
        );

    return true;
}


// ************************************************************************* //
