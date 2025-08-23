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
    (c) 2023-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "solidBodyMotionFunctions/crossWind/crossWind.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(crossWind, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        crossWind,
        registry
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::crossWind::crossWind
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    solidBodyMotionFunction(obr, SBMFCoeffs, frameName),
    velocity_(nullptr),
    frameName_(SBMFCoeffs.lookup<word>("frameForRotation")),
    frameForRotation_
    (
        coordinateFrame::New(dynamic_cast<const fvMesh&>(obr), frameName_)
    )
{
    crossWind::read(SBMFCoeffs);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::crossWind::~crossWind()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion Foam::solidBodyMotionFunctions::crossWind::transformation
(
    const scalar t1,
    const scalar t2
) const
{
    frameForRotation_.updateState();

    septernion TR(septernion::I);

    forAll(frameForRotation_.parents(), i)
    {
        if (!frameForRotation_.parents()[i].isDynamic())
        {
            TR *=
                inv(frameForRotation_.parents()[i].decoupledTransformation());
        }
    }

    DebugInFunction
        << "Transformation from time " << t1
        << " to " << t2 << " is " << TR << endl;

    return TR;
}


Foam::vectorTuple Foam::solidBodyMotionFunctions::crossWind::velocity() const
{
    vector vel(velocity_->value(time_.value()));

    // Change direction of velocity vector by nested frames
    return vectorTuple(frame().cartesianSys().transform(vel), Zero);
}


bool Foam::solidBodyMotionFunctions::crossWind::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);
    velocity_.reset(Function1<vector>::New("velocity", SBMFCoeffs_).ptr());

    return true;
}


// ************************************************************************* //
