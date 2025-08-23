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

#include "solidBodyMotionFunctions/tabulated6DoFMotion/tabulated6DoFMotion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "primitives/Tuple2/Tuple2.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "interpolations/interpolateSplineXY/interpolateSplineXY.H"
#include "global/constants/mathematical/mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(tabulated6DoFMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        tabulated6DoFMotion,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        tabulated6DoFMotion,
        registry
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::tabulated6DoFMotion::tabulated6DoFMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime,
    const word& frameName
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime, frameName),
    origin_(SBMFCoeffs_.lookup("origin"))
{
    tabulated6DoFMotion::read(SBMFCoeffs);
}

Foam::solidBodyMotionFunctions::tabulated6DoFMotion::tabulated6DoFMotion
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    solidBodyMotionFunction(obr, SBMFCoeffs, frameName),
    origin_(CofR0())
{
    tabulated6DoFMotion::read(SBMFCoeffs);
    setIncrementalMotion(tabulated6DoFMotion::typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::tabulated6DoFMotion::~tabulated6DoFMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::tabulated6DoFMotion::transformation
(
    const scalar t1,
    const scalar t2
) const
{
    if (t2 < times_[0])
    {
        FatalErrorInFunction
            << "current time (" << t2
            << ") is less than the minimum in the data table ("
            << times_[0] << ')'
            << exit(FatalError);
    }

    if (t2 > times_.last())
    {
        FatalErrorInFunction
            << "current time (" << t2
            << ") is greater than the maximum in the data table ("
            << times_.last() << ')'
            << exit(FatalError);
    }

    vector2DVector TRV = interpolateSplineXY(t2, times_, values_);

    // If the frame is available motion is only incremental
    if (hasFrame())
    {
        TRV -= interpolateSplineXY(t1, times_, values_);

        TRV[0] = frame().cartesianSys0().transform(TRV[0]);
        TRV[1] = frame().cartesianSys0().transform(TRV[1]);

        // Update origin and axis if relative motion
        origin_ = frame().coorSys().origin();
    }

    quaternion R(quaternion::XYZ, TRV[1]*pi/180.0);
    septernion TR(septernion(-origin_ - TRV[0])*R*septernion(origin_));

    DebugInFunction
        << "Transformation from time " << t1
        << " to " << t2 << " is " << TR << endl;

    return TR;
}


Foam::vectorTuple
Foam::solidBodyMotionFunctions::tabulated6DoFMotion::velocity() const
{
    const scalar t = time_.value();
    const scalar dt = time_.deltaTValue();

    vector2DVector TRV =
        interpolateSplineXY(t, times_, values_)
      - interpolateSplineXY(t - dt, times_, values_);

    if (hasFrame())
    {
        TRV[0] = frame().cartesianSys().transform(TRV[0]);
        TRV[1] = frame().cartesianSys().transform(TRV[1]);
    }

    return vectorTuple(TRV[0]/dt, (TRV[1]*pi/180.0)/dt);
}


bool Foam::solidBodyMotionFunctions::tabulated6DoFMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    // This mess has to be replaced by function1 as it was done in,
    // however this still require quite a bit of re-structure of Function1
    fileName newTimeDataFileName
    (
        fileName(SBMFCoeffs_.lookup("timeDataFileName")).expand()
    );

    if (newTimeDataFileName != timeDataFileName_)
    {
        timeDataFileName_ = newTimeDataFileName;
        IFstream dataStream(timeDataFileName_);

        if (dataStream.good())
        {
            List<Tuple2<scalar, vector2DVector>> timeValues(dataStream);

            times_.setSize(timeValues.size());
            values_.setSize(timeValues.size());

            forAll(timeValues, i)
            {
                times_[i] = timeValues[i].first();
                values_[i] = timeValues[i].second();
            }
        }
        else
        {
            FatalErrorInFunction
                << "Cannot open time data file " << timeDataFileName_
                << exit(FatalError);
        }
    }

    return true;
}


// ************************************************************************* //
