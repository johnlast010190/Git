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
    (c) 2019-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "solidBodyMotionFunctions/oscillatingRotatingMotion/oscillatingRotatingMotion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "global/unitConversion/unitConversion.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(oscillatingRotatingMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        oscillatingRotatingMotion,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        oscillatingRotatingMotion,
        registry
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::oscillatingRotatingMotion::
oscillatingRotatingMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime,
    const word& frameName
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime, frameName)
{
    oscillatingRotatingMotion::read(SBMFCoeffs);
}


Foam::solidBodyMotionFunctions::oscillatingRotatingMotion::
oscillatingRotatingMotion
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    solidBodyMotionFunction(obr, SBMFCoeffs, frameName),
    origin_(CofR0())
{
    oscillatingRotatingMotion::read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::oscillatingRotatingMotion::
~oscillatingRotatingMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::oscillatingRotatingMotion::transformation
(
    const scalar t1,
    const scalar t2
) const
{
    vector eulerAngles = amplitude_->value(t2)*sin(omega_->integrate(0.0, t2));

    // Update origin and axis if relative motion
    if (hasFrame())
    {
        if (isIncrementalMotion())
        {
            eulerAngles -=
                amplitude_->value(t1)*sin(omega_->integrate(0.0, t1));
            eulerAngles = frame().cartesianSys().transform(eulerAngles);
            origin_ = frame().coorSys().origin();
        }
        else
        {
            eulerAngles = frame().cartesianSys0().transform(eulerAngles);
        }
    }

    // Convert the rotational motion from deg to rad
    eulerAngles *= pi/180.0;

    const quaternion R(quaternion::XYZ, eulerAngles);
    const septernion TR(septernion(-origin_)*R*septernion(origin_));

    DebugInFunction
        << "Transformation from time " << t1
        << " to " << t2 << " is " << TR << endl;

    return TR;
}



Foam::vectorTuple
Foam::solidBodyMotionFunctions::oscillatingRotatingMotion::velocity() const
{
    const scalar t = time_.value();
    const scalar intOmegaT = omega_->integrate(0, t);
    vector omega
    (
        amplitude_->derivative(t)*sin(intOmegaT)
      + amplitude_->value(t)*omega_->value(t)*cos(intOmegaT)
    );
    if (hasFrame())
    {
        omega = frame().cartesianSys().transform(omega);
    }

    return vectorTuple(Zero, (omega*pi/180.0));
}


bool Foam::solidBodyMotionFunctions::oscillatingRotatingMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    if (SBMFCoeffs_.found("origin"))
    {
        origin_ = SBMFCoeffs_.lookup<vector>("origin");
    }
    omega_.reset(Function1<scalar>::New("omega", SBMFCoeffs_));
    amplitude_.reset(Function1<vector>::New("amplitude", SBMFCoeffs_));

    return true;
}


// ************************************************************************* //
