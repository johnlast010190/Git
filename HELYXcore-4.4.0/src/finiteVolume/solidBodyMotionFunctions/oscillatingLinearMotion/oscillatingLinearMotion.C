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

#include "solidBodyMotionFunctions/oscillatingLinearMotion/oscillatingLinearMotion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(oscillatingLinearMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        oscillatingLinearMotion,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        oscillatingLinearMotion,
        registry
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::oscillatingLinearMotion::oscillatingLinearMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime,
    const word& frameName
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime, frameName)
{
    oscillatingLinearMotion::read(SBMFCoeffs);
}

Foam::solidBodyMotionFunctions::oscillatingLinearMotion::oscillatingLinearMotion
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    solidBodyMotionFunction(obr, SBMFCoeffs, frameName)
{
    oscillatingLinearMotion::read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::oscillatingLinearMotion::
~oscillatingLinearMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::oscillatingLinearMotion::transformation
(
    const scalar t1,
    const scalar t2
) const
{
    vector displacement =
        amplitude_->value(t2)*sin(omega_->integrate(0.0, t2));

    // Update origin and axis if relative motion
    if (hasFrame())
    {
        if (isIncrementalMotion())
        {
            displacement -=
                amplitude_->value(t1)*sin(omega_->integrate(0.0, t1));
            displacement = frame().cartesianSys().transform(displacement);
        }
        else
        {
            displacement = frame().cartesianSys0().transform(displacement);
        }
    }

    quaternion R(1);
    septernion TR(septernion(-displacement)*R);

    DebugInFunction
        << "Transformation from time " << t1
        << " to " << t2 << " is " << TR << endl;

    return TR;
}


Foam::vectorTuple
Foam::solidBodyMotionFunctions::oscillatingLinearMotion::velocity() const
{
    const scalar t = time_.value();
    const scalar intOmegaT = omega_->integrate(0, t);
    vector velocity
    (
        amplitude_->derivative(t)*sin(intOmegaT)
      + amplitude_->value(t)*omega_->value(t)*cos(intOmegaT)
    );
    if (hasFrame())
    {
        velocity = frame().cartesianSys().transform(velocity);
    }

    return vectorTuple(velocity, Zero);
}


bool Foam::solidBodyMotionFunctions::oscillatingLinearMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);
    omega_.reset(Function1<scalar>::New("omega", SBMFCoeffs_));
    amplitude_.reset(Function1<vector>::New("amplitude", SBMFCoeffs_));

    return true;
}


// ************************************************************************* //
