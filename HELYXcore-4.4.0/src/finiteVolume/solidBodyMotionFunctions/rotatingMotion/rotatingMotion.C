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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2019-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "solidBodyMotionFunctions/rotatingMotion/rotatingMotion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(rotatingMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        rotatingMotion,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        rotatingMotion,
        registry
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::rotatingMotion::rotatingMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime,
    const word& frameName
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime, frameName),
    origin_(SBMFCoeffs_.lookup("origin")),
    axis_(SBMFCoeffs_.lookup("axis")),
    omega_
    (
        SBMFCoeffs_.found("omega")
      ? Function1<scalar>::New("omega", SBMFCoeffs_)
      : nullptr
    ),
    angle_
    (
        !SBMFCoeffs_.found("omega")
       ? Function1<scalar>::New("angle", SBMFCoeffs_)
       : nullptr
    )
{}


Foam::solidBodyMotionFunctions::rotatingMotion::rotatingMotion
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    solidBodyMotionFunction(obr, SBMFCoeffs, frameName),
    origin_(CofR0()),
    axis_(normalised(frame().axis0())),
    omega_
    (
        SBMFCoeffs_.found("omega")
      ? Function1<scalar>::New("omega", SBMFCoeffs_)
      : nullptr
    ),
    angle_
    (
        !SBMFCoeffs_.found("omega")
       ? Function1<scalar>::New("angle", SBMFCoeffs_)
       : nullptr
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::rotatingMotion::~rotatingMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::rotatingMotion::transformation
(
    const scalar t1,
    const scalar t2
) const
{
    // Update origin and axis if relative motion
    if (hasFrame() && isIncrementalMotion())
    {
        origin_ = frame().coorSys().origin();
        axis_ = normalised(frame().axis());
    }

    // Rotation around axis
    const scalar angle =
        omega_.valid()
      ? omega_->integrate(t1, t2)
      : (angle_->value(t2) - angle_->value(t1));
    const quaternion R(axis_, angle);
    const septernion TR(septernion(-origin_)*R*septernion(origin_));

    DebugInFunction
        << "Time = " << t2 << ", t1 = " << t1
        << ", transformation: " << TR << endl;

    return TR;
}


Foam::vectorTuple
Foam::solidBodyMotionFunctions::rotatingMotion::velocity() const
{
    const scalar time = time_.value();
    const scalar deltaT = time_.deltaT().value();
    const scalar prevTime = time > 0 ? (time - deltaT) : 0;

    // Note: Derivative might be better but it might not work for "REPEAT"
    return
        vectorTuple
        (
            Zero,
            normalised(hasFrame() ? frame().axis() : axis_)
           *(
                omega_.valid()
              ? omega_->value(time)
              : ((angle_->value(time) - angle_->value(prevTime))/deltaT)
            )
        );
}


bool Foam::solidBodyMotionFunctions::rotatingMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);
    if (SBMFCoeffs_.found("omega"))
    {
        omega_.reset(Function1<scalar>::New("omega", SBMFCoeffs_).ptr());
        angle_.reset(nullptr);
    }
    else
    {
        omega_.reset(nullptr);
        angle_.reset(Function1<scalar>::New("angle", SBMFCoeffs_).ptr());
    }
    return true;
}


// ************************************************************************* //
