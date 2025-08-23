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

\*---------------------------------------------------------------------------*/

#include "exponential/exponential.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{
    defineTypeNameAndDebug(exponential, 0);
    addToRunTimeSelectionTable(distributionModel, exponential, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionModels::exponential::exponential
(
    const dictionary& dict,
    Random& rndGen
)
:
    distributionModel(typeName, dict, rndGen),
    minValue_(distributionModelDict_.lookup<scalar>("minValue")),
    maxValue_(distributionModelDict_.lookup<scalar>("maxValue")),
    lambda_(distributionModelDict_.lookup<scalar>("lambda"))
{
    check();
}


Foam::distributionModels::exponential::exponential(const exponential& p)
:
    distributionModel(p),
    minValue_(p.minValue_),
    maxValue_(p.maxValue_),
    lambda_(p.lambda_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributionModels::exponential::~exponential()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::exponential::sample() const
{
    scalar y = rndGen_.sample01<scalar>();
    scalar K = exp(-lambda_*maxValue_) - exp(-lambda_*minValue_);
    return -(1.0/lambda_)*log(exp(-lambda_*minValue_) + y*K);
}


Foam::scalar Foam::distributionModels::exponential::minValue() const
{
    return minValue_;
}


Foam::scalar Foam::distributionModels::exponential::maxValue() const
{
    return maxValue_;
}


Foam::scalar Foam::distributionModels::exponential::meanValue() const
{
    return 1.0/lambda_;
}


// ************************************************************************* //
