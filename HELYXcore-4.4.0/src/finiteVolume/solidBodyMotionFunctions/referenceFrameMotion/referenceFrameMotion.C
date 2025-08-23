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
#include "solidBodyMotionFunctions/referenceFrameMotion/referenceFrameMotion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(referenceFrameMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        referenceFrameMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::referenceFrameMotion::referenceFrameMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime,
    const word& frameName
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime, frameName)
{
    DeprecationWarningInFunction
    (
        "referenceFrameMotion",
        "motion type",
        40200,
        "Please the frame definition directly in dynamicMeshDict."
    );

    referenceFrameMotion::read(SBMFCoeffs);
    // this motion function is always dynamic
    coorFramePtr_->resetDynamic(true);

    // This whole class should be removed ones it is removed from GUI
    if (coorFramePtr_->isIncrementalMotion())
    {
        incrementalMotion_ = true;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::referenceFrameMotion::~referenceFrameMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::referenceFrameMotion::transformation() const
{
    coorFramePtr_->updateState();

    return coorFramePtr_->transformation();
}


Foam::septernion
Foam::solidBodyMotionFunctions::referenceFrameMotion::transformation
(
    const scalar t1,
    const scalar t2
) const
{
    NotImplemented;
}


Foam::vectorTuple
Foam::solidBodyMotionFunctions::referenceFrameMotion::velocity() const
{
    coorFramePtr_->updateState();
    return coorFramePtr_->velocity();
}


bool Foam::solidBodyMotionFunctions::referenceFrameMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    frameName_ = SBMFCoeffs_.lookup<word>("referenceFrame");

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>
        (
            SBMFCoeffs_.lookupOrDefault<word>("region", "region0")
        );

    if (mesh.foundObject<coordinateFrame>(frameName_))
    {
        coorFramePtr_ =
            &mesh.lookupObjectRef<coordinateFrame>(frameName_);
    }
    else if (frameName_ != word::null)
    {
        coorFramePtr_ = &coordinateFrame::New(mesh, frameName_);
    }

    return true;
}


// ************************************************************************* //
