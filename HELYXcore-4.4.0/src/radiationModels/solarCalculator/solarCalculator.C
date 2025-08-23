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
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "solarCalculator/solarCalculator.H"
#include "db/Time/Time.H"
#include "global/unitConversion/unitConversion.H"
#include "global/constants/constants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solarCalculator, 0);

    template<>
    const char* NamedEnum
    <
        solarCalculator::sunDirModel,
        2
    >::names[] =
    {
        "sunDirConstant",
        "sunDirTracking"
    };

    template<>
    const char* NamedEnum
    <
        solarCalculator::sunLModel,
        4
    >::names[] =
    {
        "sunLoadConstant",
        "sunLoadTimeDependent",
        "sunLoadFairWeatherConditions",
        "sunLoadTheoreticalMaximum"
    };
}

const Foam::NamedEnum<Foam::solarCalculator::sunDirModel, 2>
  Foam::solarCalculator::sunDirectionModelTypeNames_;

const Foam::NamedEnum<Foam::solarCalculator::sunLModel, 4>
   Foam::solarCalculator::sunLoadModelTypeNames_;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solarCalculator::calculateBetaTheta()
{
    scalar runTime = 0.0;
    switch (sunDirectionModel_)
    {
        case mSunDirTracking:
        {
            runTime = mesh_.time().value();
            break;
        }
        case mSunDirConstant:
        {
            break;
        }
    }

    scalar LSM = 15.0*(dict_.lookup<scalar>("localStandardMeridian"));

    scalar D = dict_.lookup<scalar>("startDay") + runTime/86400.0;
    scalar M = 6.24004 + 0.0172*D;
    scalar EOT = -7.659*sin(M) + 9.863*sin(2*M + 3.5932);

    startTime_ = dict_.lookup<scalar>("startTime");
    scalar LST =  startTime_ + runTime/3600.0;

    scalar LON = dict_.lookup<scalar>("longitude");

    scalar AST = LST + EOT/60.0 + (LON - LSM)/15;

    scalar delta = 23.45*sin(degToRad((360*(284 + D))/365));

    scalar H = degToRad(15*(AST - 12));

    scalar L = degToRad(dict_.lookup<scalar>("latitude"));

    scalar deltaRad = degToRad(delta);
    beta_ = max(asin(cos(L)*cos(deltaRad)*cos(H) + sin(L)*sin(deltaRad)), 1e-3);
    theta_ = acos((sin(beta_)*sin(L) - sin(deltaRad))/(cos(beta_)*cos(L)));

    // theta is the angle between the SOUTH axis and the Sun
    // If the hour angle is lower than zero (morning) the Sun is positioned
    // on the East side.
    if (H < 0)
    {
        theta_ += 2*(constant::mathematical::pi - theta_);
    }

    if (debug)
    {
        Info<< tab << "altitude : " << radToDeg(beta_) << endl;
        Info<< tab << "azimuth  : " << radToDeg(theta_) << endl;
    }
}


void Foam::solarCalculator::calculateSunDirection()
{
    coord_.reset
    (
        new coordinateSystem("grid", Zero, gridUp_, eastDir_)
    );

    // Assuming 'z' vertical, 'y' North and 'x' East
    direction_.z() = -sin(beta_);
    direction_.y() =  cos(beta_)*cos(theta_); // South axis
    direction_.x() =  cos(beta_)*sin(theta_); // West axis

    direction_ /= mag(direction_);

    if (debug)
    {
        Info<< "Sun direction in absolute coordinates : " << direction_ <<endl;
    }

    // Transform to actual coordinate system
    direction_ = coord_->transform(direction_);

    if (debug)
    {
        Info<< "Sun direction in the Grid coordinates : " << direction_ <<endl;
    }
}


void Foam::solarCalculator::init()
{
    switch (sunDirectionModel_)
    {
        case mSunDirConstant:
        {
            if (dict_.found("sunDirection"))
            {
                direction_ = dict_.lookup<vector>("sunDirection");
                direction_ /= mag(direction_);
            }
            else
            {
                updateOrientation();
                calculateBetaTheta();
                calculateSunDirection();
            }

            break;
        }
        case mSunDirTracking:
        {
            if (word(mesh_.schemes().ddtScheme("default")) == "steadyState")
            {
                FatalErrorInFunction
                    << " Sun direction model can not be sunDirtracking if the "
                    << " case is steady " << nl << exit(FatalError);
            }

            sunTrackingUpdateInterval_ =
                dict_.lookup<scalar>("sunTrackingUpdateInterval");

            updateOrientation();
            calculateBetaTheta();
            calculateSunDirection();
            break;
        }
    }

    switch (sunLoadModel_)
    {
        case mSunLoadConstant:
        {
            directSolarRad_ = dict_.lookup<scalar>("directSolarRad");
            diffuseSolarRad_ = dict_.lookup<scalar>("diffuseSolarRad");
            break;
        }
        case mSunLoadTimeDependent:
        {
            directSolarRadPtr_.reset
            (
                Function1<scalar>::New
                (
                    "directSolarRad",
                    dict_
                ).ptr()
            );
            directSolarRad_ =
            directSolarRadPtr_->value(mesh_.time().timeOutputValue());
            diffuseSolarRadPtr_.reset
            (
                Function1<scalar>::New
                (
                    "diffuseSolarRad",
                    dict_
                ).ptr()
            );
            diffuseSolarRad_ =
            diffuseSolarRadPtr_->value(mesh_.time().timeOutputValue());
            break;
        }
        case mSunLoadFairWeatherConditions:
        {
            dict_.readIfPresent
            (
                "skyCloudCoverFraction",
                skyCloudCoverFraction_
            );

            A_ = dict_.lookup<scalar>("A");
            B_ = dict_.lookup<scalar>("B");

            if (dict_.found("beta"))
            {
                beta_ = dict_.lookup<scalar>("beta");
            }
            else
            {
                calculateBetaTheta();
            }

            directSolarRad_ =
                (1.0 - 0.75*pow(skyCloudCoverFraction_, 3.0))
              * A_/exp(B_/sin(beta_));

            groundReflectivity_ =
                dict_.lookup<scalar>("groundReflectivity");

            break;
        }
        case mSunLoadTheoreticalMaximum:
        {
            Setrn_ = dict_.lookup<scalar>("Setrn");
            SunPrime_ = dict_.lookup<scalar>("SunPrime");
            directSolarRad_ = Setrn_*SunPrime_;

            groundReflectivity_ =
                dict_.lookup<scalar>("groundReflectivity");
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solarCalculator::solarCalculator
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    dict_(dict),
    direction_(Zero),
    directSolarRad_(0.0),
    directSolarRadPtr_(),
    diffuseSolarRad_(0.0),
    diffuseSolarRadPtr_(),
    groundReflectivity_(0.0),
    A_(0.0),
    B_(0.0),
    beta_(0.0),
    theta_(0.0),
    skyCloudCoverFraction_(0.0),
    Setrn_(0.0),
    SunPrime_(0.0),
    C_(dict.lookupOrDefault<scalar>("C", 0.058)),
    sunDirectionModel_
    (
        sunDirectionModelTypeNames_.read(dict.lookup("sunDirectionModel"))
    ),
    sunLoadModel_
    (
        sunLoadModelTypeNames_.read(dict.lookup("sunLoadModel"))
    ),
    coord_()
{
    init();
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solarCalculator::~solarCalculator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solarCalculator::correctSunDirection()
{
    switch (sunDirectionModel_)
    {
        case mSunDirConstant:
        {
            break;
        }
        case mSunDirTracking:
        {
            calculateBetaTheta();
            calculateSunDirection();
            directSolarRad_ = A_/exp(B_/sin(max(beta_, ROOTVSMALL)));
            break;
        }
    }
}


void Foam::solarCalculator::correctSunFlux()
{
    if (sunLoadModel_==mSunLoadTimeDependent)
    {
        if (!diffuseSolarRadPtr_)
        {
            FatalErrorInFunction<< "DiffuseSolarRad does not exist."
                                << exit(FatalError);
        }
        diffuseSolarRad_ =
        diffuseSolarRadPtr_->value(mesh_.time().timeOutputValue());

        if (!directSolarRadPtr_)
        {
            FatalErrorInFunction<< "DirectSolarRad does not exist."
                                << exit(FatalError);
        }
        directSolarRad_ =
        directSolarRadPtr_->value(mesh_.time().timeOutputValue());
    }
}


void Foam::solarCalculator::updateOrientation()
{
    gridUp_ = dict_.lookup<vector>("gridUp");
    gridUp_ /= mag(gridUp_);

    eastDir_ = dict_.lookup<vector>("gridEast");
    eastDir_ /= mag(eastDir_);
}

// ************************************************************************* //
