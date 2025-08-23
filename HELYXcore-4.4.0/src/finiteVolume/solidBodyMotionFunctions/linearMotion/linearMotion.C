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
    (c) 2020-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "solidBodyMotionFunctions/linearMotion/linearMotion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(linearMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        linearMotion,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        linearMotion,
        registry
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::linearMotion::linearMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime,
    const word& frameName
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime, frameName),
    velocity_(nullptr)
{
    linearMotion::read(SBMFCoeffs);
}


Foam::solidBodyMotionFunctions::linearMotion::linearMotion
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    solidBodyMotionFunction(obr, SBMFCoeffs, frameName),
    velocity_(nullptr)
{
    linearMotion::read(SBMFCoeffs);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::linearMotion::~linearMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::linearMotion::transformation
(
    const scalar t1,
    const scalar t2
) const
{
    // Translation of centre of gravity with constant velocity
    vector displacement = velocity_->integrate(t1, t2);

    if (hasFrame())
    {
        if (isIncrementalMotion())
        {
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
Foam::solidBodyMotionFunctions::linearMotion::velocity() const
{
    vector vel(velocity_->value(time_.value()));

    if (hasFrame())
    {
        vel = frame().cartesianSys().transform(vel);
    }

    return vectorTuple(vel, Zero);
}


bool Foam::solidBodyMotionFunctions::linearMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);
    velocity_.reset(Function1<vector>::New("velocity", SBMFCoeffs_).ptr());

    return true;
}


// ************************************************************************* //
