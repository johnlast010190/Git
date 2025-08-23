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
    (c) 2023-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "solidBodyMotionFunctions/PIDSpeedRegulatorMotion/PIDSpeedRegulatorMotion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "interpolations/interpolateXY/interpolateXY.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(PIDSpeedRegulatorMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        PIDSpeedRegulatorMotion,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        PIDSpeedRegulatorMotion,
        registry
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::PIDSpeedRegulatorMotion::PIDSpeedRegulatorMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime,
    const word& frameName
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime, frameName),
    angle_(0.0),
    origin_(SBMFCoeffs_.lookup("origin")),
    axis_(SBMFCoeffs_.lookup("axis"))
{
    PIDSpeedRegulatorMotion::read(SBMFCoeffs);
    setIncrementalMotion(PIDSpeedRegulatorMotion::typeName);
}


Foam::solidBodyMotionFunctions::PIDSpeedRegulatorMotion::PIDSpeedRegulatorMotion
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    solidBodyMotionFunction(obr, SBMFCoeffs, frameName),
    angle_(0.0),
    origin_(coorFramePtr_->coorSys().origin()),
    axis_(coorFramePtr_->axis())
{
    PIDSpeedRegulatorMotion::read(SBMFCoeffs);
    setIncrementalMotion(PIDSpeedRegulatorMotion::typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::PIDSpeedRegulatorMotion::
~PIDSpeedRegulatorMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solidBodyMotionFunctions::PIDSpeedRegulatorMotion::
updateOriginAndAxis() const
{
    // Update origin and axis if relative motion
    if (hasFrame() && isIncrementalMotion())
    {
        origin_ = frame().coorSys().origin();
        axis_ = normalised(frame().axis());
    }
}


Foam::scalar
Foam::solidBodyMotionFunctions::PIDSpeedRegulatorMotion::torque() const
{
    const functionObjectList& functions = time_.functionObjects();
    const IOdictionary& stateDict = functions.stateDict();

    vector normalMoment = vector::zero;
    vector tangentialMoment = vector::zero;

    if (stateDict.found("results"))
    {
        const dictionary& resultsDict = stateDict.subDict("results");

        if (resultsDict.found(functionName_))
        {
            const dictionary& objectDict = resultsDict.subDict(functionName_);

            if (objectDict.found("vector"))
            {
                const dictionary& resultTypeDict =
                    objectDict.subDict("vector");

                resultTypeDict.readIfPresent<vector>
                (
                    "normalMoment",
                    normalMoment
                );

                resultTypeDict.readIfPresent<vector>
                (
                    "tangentialMoment",
                    tangentialMoment
                );
            }
        }
    }

    vector totalMoment = normalMoment + tangentialMoment;
    updateOriginAndAxis();
    return (totalMoment & axis_);
}


Foam::septernion
Foam::solidBodyMotionFunctions::PIDSpeedRegulatorMotion::transformation
(
    const scalar t1,
    const scalar t2
) const
{
    scalar dt = t2 - t1;
    scalar currentTorque = torque();

    const regulatorItem& last = results_.last();
    const scalar currentOmega = last.regulatedSpeed();
    const scalar currentLimit = limitCurve_->value(currentOmega);

    const label nItems = results_.size();

    if (nItems > 2 && t2 >= startTime_)
    {
        // Set propellor curve
        scalarField xVal(12);
        scalarField yVal(12);

        for (label i = 0; i < 12; i++)
        {
            scalar d = i * scalar(0.1);

            xVal[i] = currentOmega * d;
            yVal[i] = currentTorque*Foam::pow(xVal[i]/currentOmega, 3);
        }

        const scalar targetOmega =
            interpolateXY
            (
                (
                    currentTorque > currentLimit
                  ? currentLimit
                  : designPoint_.second()
                ),
                yVal,
                xVal
            );

        scalarField Is(3, Zero);
        scalarField Should(3, Zero);
        scalarField Differences(3, Zero);

        label resultI = 0;
        forAllReverse(results_, i)
        {
            const regulatorItem& item = results_[i];
            if (resultI == 0)
            {
                Is[resultI] = currentOmega;
                Should[resultI] = targetOmega;
            }
            else
            {
                Is[resultI] = item.propellerSpeed();
                Should[resultI] = item.regulatedSpeed();
            }
            Differences[resultI] = Should[resultI] - Is[resultI];

            resultI++;
            if (resultI == 3)
            {
                break;
            }
        }

        // Actual PID formula, calculating new propeller omega
        omega_ =
            Is[2]
          + (
                kP_*Differences[2]
              + tV_*(Differences[2] - Differences[1])
              + tN_*sum(Differences)
            );
    }

    omega_ = omega_ > speedLimits_.second() ? speedLimits_.second() : omega_;
    omega_ = omega_ < speedLimits_.first() ? speedLimits_.first() : omega_;

    const regulatorItem& lastItem = results_.last();

    scalar newSteering =
        lastItem.steeringAngle()
      + (steeringSpeed_/(2*constant::mathematical::pi))*dt;

    Info<< "PID Speed Regulator : " << nl
        << "Steering : " << newSteering
        << " Current Omega : " << currentOmega << nl
        << "Current Torque : " << currentTorque
        << " Current Limit : " << currentLimit
        << " Regulated  Speed : " << omega_ << nl
        << endl;

    results_.append
    (
        regulatorItem
        (
            newSteering,
            currentOmega,
            currentTorque,
            currentLimit,
            omega_
        )
    );

    // Convert from revolutions to radians and increment the angle
    angle_ += dt*omega_;

    quaternion R(axis_, angle_);
    septernion TR(septernion(-origin_)*R*septernion(origin_));

    DebugInFunction
        << "Transformation from time " << t1
        << " to " << t2 << " is " << TR << endl;

    return TR;
}


Foam::vectorTuple
Foam::solidBodyMotionFunctions::PIDSpeedRegulatorMotion::velocity() const
{
    const scalar t  = time_.value();
    scalar currentTorque = torque();

    const regulatorItem& last = results_.last();
    const scalar currentOmega = last.regulatedSpeed();
    const scalar currentLimit = limitCurve_->value(currentOmega);

    scalar omega(Zero);

    label nItems = results_.size();

    if (nItems > 2 && t >= startTime_)
    {
        // Set propellor curve
        scalarField xVal(12);
        scalarField yVal(12);

        for (label i = 0; i < 12; i++)
        {
            const scalar d = scalar(i)*0.1;
            xVal[i] = currentOmega*d;
            yVal[i] = currentTorque*Foam::pow(xVal[i]/currentOmega, 3);
        }

        const scalar targetOmega =
            interpolateXY
            (
                (
                    currentTorque > currentLimit
                  ? currentLimit
                  : designPoint_.second()
                ),
                yVal,
                xVal
            );

        scalarField Is(3, Zero);
        scalarField Should(3, Zero);
        scalarField Differences(3, Zero);

        label resultI = 0;
        forAllReverse(results_, i)
        {
            const regulatorItem& item = results_[i];
            if (resultI == 0)
            {
                Is[resultI] = currentOmega;
                Should[resultI] = targetOmega;
            }
            else
            {
                Is[resultI] = item.propellerSpeed();
                Should[resultI] = item.regulatedSpeed();
            }
            Differences[resultI] = Should[resultI] - Is[resultI];

            resultI++;
            if (resultI == 3)
            {
                break;
            }
        }

        // Actual PID formula, calculating new propeller omega
        omega =
            Is[2]
          + (
                kP_*Differences[2]
              + tV_*(Differences[2] - Differences[1])
              + tN_ * sum(Differences)
            );
    }

    omega = omega > speedLimits_.second() ? speedLimits_.second() : omega;
    omega = omega < speedLimits_.first() ? speedLimits_.first() : omega;

    return vectorTuple(Zero, axis_*omega);
}


bool Foam::solidBodyMotionFunctions::PIDSpeedRegulatorMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    if (mag(axis_) > SMALL)
    {
        axis_ /= mag(axis_);
    }
    else
    {
        FatalErrorInFunction
            << "Keyword axis set incorrectly: " << axis_ << endl;
    }

    // Forces function object name for calculation of torque
    functionName_ = SBMFCoeffs_.lookup<word>("functionName");

    // Lower and upper propeller speed limits
    speedLimits_ = SBMFCoeffs_.lookup<Pair<scalar>>("speedLimits");

    // Regulator parameters
    kP_ = SBMFCoeffs_.lookup<scalar>("KP");
    tN_ = SBMFCoeffs_.lookup<scalar>("TN");
    tV_ = SBMFCoeffs_.lookup<scalar>("TV");
    steeringSpeed_ = SBMFCoeffs_.lookup<scalar>("steeringSpeed");
    designPoint_ = SBMFCoeffs_.lookup<Pair<scalar>>("designPoint");

    // Set initial Omega to design point
    omega_ = designPoint_.first();

    limitCurve_ = Function1<scalar>::New("limitCurve", SBMFCoeffs_);

    startTime_ =
        SBMFCoeffs_.lookupOrDefault<scalar>("regulateStartTime", 0.0);

    results_.append
    (
        regulatorItem
        (
            scalar(0),
            designPoint_.first(), // Design omega
            designPoint_.second(), // Design torque
            limitCurve_->value(designPoint_.first()),
            designPoint_.first()
         )
    );

    return true;
}


// ************************************************************************* //
