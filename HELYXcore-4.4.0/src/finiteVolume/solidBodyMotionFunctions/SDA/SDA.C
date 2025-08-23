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

#include "solidBodyMotionFunctions/SDA/SDA.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/constants/mathematical/mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(SDA, 0);
    addToRunTimeSelectionTable(solidBodyMotionFunction, SDA, dictionary);
    addToRunTimeSelectionTable(solidBodyMotionFunction, SDA, registry);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::SDA::SDA
(
    const dictionary& SBMFCoeffs,
    const Time& runTime,
    const word& frameName
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime, frameName),
    origin_(SBMFCoeffs_.lookup("CofG"))
{
    SDA::read(SBMFCoeffs);
}


Foam::solidBodyMotionFunctions::SDA::SDA
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    solidBodyMotionFunction(obr, SBMFCoeffs, frameName),
    origin_
    (
        SBMFCoeffs_.lookupOrAddDefault<vector>
        (
            "CofG",
            frame().coorSys0().origin()
        )
    )
{
    SDA::read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::SDA::~SDA()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Vector2D<Foam::vector>
Foam::solidBodyMotionFunctions::SDA::transformationAndRotation
(
    const scalar t
) const
{
    // Current roll period [sec]
    const scalar Tpi = Tp_ + dTp_*(t/dTi_);

    // Current Freq [/sec]
    const scalar wr = twoPi/Tpi;

    // Current Phase for roll [rad]
    const scalar r = dTp_/dTi_;
    const scalar u = Tp_ + r*t;
    const scalar phr = twoPi*((Tp_/u - 1) + log(mag(u)) - log(Tp_))/r;

    // Current Phase for Sway [rad]
    const scalar phs = phr + pi;

    // Current Phase for Heave [rad]
    const scalar phh = phr + piByTwo;

    const scalar rollA =
        max(rollAmax_*exp(-sqr(Tpi - Tpn_)/(2*Q_)), rollAmin_);

    const vector T
    (
        0,
        swayA_*(sin(wr*t + phs) - sin(phs)),
        heaveA_*(sin(wr*t + phh) - sin(phh))
    );

    return Vector2D<vector>(T, vector(rollA*sin(wr*t + phr), 0, 0));
}


Foam::septernion Foam::solidBodyMotionFunctions::SDA::transformation
(
    const scalar t1,
    const scalar t2
) const
{
    Vector2D<vector> TRV = transformationAndRotation(t2);
    if (hasFrame())
    {
        if (isIncrementalMotion())
        {
            TRV -= transformationAndRotation(t1);
        }

        // Motion is always defined in initial frame => probably doesn't
        // make sense to define in reference to moving frame
        TRV[0] = frame().cartesianSys0().transform(TRV[0]);
        TRV[1] = frame().cartesianSys0().transform(TRV[1]);
    }
    quaternion R(quaternion::XYZ, TRV[1]);
    septernion TR(septernion(-origin_ - TRV[0])*R*septernion(origin_));

    DebugInFunction
        << "Transformation from time " << t1
        << " to " << t2 << " is " << TR << endl;

    return TR;
}


Foam::vectorTuple Foam::solidBodyMotionFunctions::SDA::velocity() const
{
    const scalar t = time_.value();
    const scalar dt = time_.deltaTValue();

    Vector2D<vector> TRV =
        transformationAndRotation(t) - transformationAndRotation(t - dt);

    if (hasFrame())
    {
        TRV[0] = frame().cartesianSys().transform(TRV[0]);
        TRV[1] = frame().cartesianSys().transform(TRV[1]);
    }

    return vectorTuple(TRV[0]/dt, (TRV[1]*pi/180.0)/dt);
}


bool Foam::solidBodyMotionFunctions::SDA::read(const dictionary& SBMFCoeffs)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    lamda_ = SBMFCoeffs_.lookup<scalar>("lamda");
    rollAmax_ = SBMFCoeffs_.lookup<scalar>("rollAmax");
    rollAmin_ = SBMFCoeffs_.lookup<scalar>("rollAmin");
    heaveA_ = SBMFCoeffs_.lookup<scalar>("heaveA");
    swayA_ = SBMFCoeffs_.lookup<scalar>("swayA");
    Q_ = SBMFCoeffs_.lookup<scalar>("Q");
    Tp_ = SBMFCoeffs_.lookup<scalar>("Tp");
    Tpn_ = SBMFCoeffs_.lookup<scalar>("Tpn");
    dTi_ = SBMFCoeffs_.lookup<scalar>("dTi");
    dTp_ = SBMFCoeffs_.lookup<scalar>("dTp");

    // Rescale parameters according to the given scale parameter
    if (lamda_ > 1 + SMALL)
    {
        heaveA_ /= lamda_;
        swayA_ /= lamda_;
        Tp_ /= sqrt(lamda_);
        Tpn_ /= sqrt(lamda_);
        dTi_ /= sqrt(lamda_);
        dTp_ /= sqrt(lamda_);
    }

    return true;
}


// ************************************************************************* //
