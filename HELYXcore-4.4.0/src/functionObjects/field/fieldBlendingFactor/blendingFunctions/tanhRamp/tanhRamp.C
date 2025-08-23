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
    (c) 2013 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "db/dictionary/dictionary.H"
#include "fieldBlendingFactor/blendingFunctions/tanhRamp/tanhRamp.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace blendingFunctions
{

defineTypeNameAndDebug(tanhRamp, 0);
addToRunTimeSelectionTable(blendingFunction, tanhRamp, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tanhRamp::tanhRamp(const dictionary& dict, bool dum)
:
    blendingFunction(dict),
    minValue_(dict.lookup<scalar>("minValue")),
    maxValue_(dict.lookup<scalar>("maxValue")),
    minCoeff_(dict.lookup<scalar>("minCoeff")),
    maxCoeff_(dict.lookup<scalar>("maxCoeff")),
    inOffset_(),
    inScale_(),
    outOffset_(),
    outScale_(),
    approxTanh_(dict.lookupOrDefault<Switch>("approxTanh", true))
{

    inOffset_ = 0.5*(minValue_ + maxValue_);
    inScale_ = 3.0/mag(maxValue_ - inOffset_);

    outOffset_ = 0.5*(minCoeff_ + maxCoeff_);
    outScale_ = maxCoeff_ - outOffset_;

    if (debug)
    {
        Info<< "Input : minValue - " << minValue_ << ", maxValue - "
             << maxValue_ << ", minCoeff - " << minCoeff_
             << ", maxCoeff - " << maxCoeff_ << endl;

        Info<< "tanh ramp coefficients: "
             << inOffset_ << ", " << inScale_ << ", "
             << outOffset_ << ", " << outScale_ <<endl;
        Info<< "Return min: " << this->operator()(minValue_) << endl;
        Info<< "Return mean: " << this->operator()(0.5*(minValue_+maxValue_))
             << endl;
        Info<< "Return max: " << this->operator()(maxValue_) << endl;
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void tanhRamp::update(scalar value)
{
    minValue_ = value;
    inOffset_ = 0.5*(minValue_ + maxValue_);
    inScale_ = 3.0/mag(maxValue_ - inOffset_);
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

scalar tanhRamp::operator()(scalar value) const
{

    scalar xin = (value - inOffset_)*inScale_;

    if (approxTanh_)
    {
        return ptanh(xin)*outScale_ + outOffset_;
    }
    else
    {
        return tanh(xin)*outScale_ + outOffset_;
    }
}


} // End namespace blendingFunctions
} // End namespace Foam

// ************************************************************************* //
