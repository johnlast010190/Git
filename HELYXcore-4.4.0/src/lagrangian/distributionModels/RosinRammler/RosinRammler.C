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

#include "RosinRammler/RosinRammler.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{
    defineTypeNameAndDebug(RosinRammler, 0);
    addToRunTimeSelectionTable(distributionModel, RosinRammler, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionModels::RosinRammler::RosinRammler
(
    const dictionary& dict,
    Random& rndGen
)
:
    distributionModel(typeName, dict, rndGen),
    minValue_(distributionModelDict_.lookup<scalar>("minValue")),
    maxValue_(distributionModelDict_.lookup<scalar>("maxValue")),
    d_(distributionModelDict_.lookup<scalar>("d")),
    n_(distributionModelDict_.lookup<scalar>("n"))
{
    check();
}


Foam::distributionModels::RosinRammler::RosinRammler(const RosinRammler& p)
:
    distributionModel(p),
    minValue_(p.minValue_),
    maxValue_(p.maxValue_),
    d_(p.d_),
    n_(p.n_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributionModels::RosinRammler::~RosinRammler()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::RosinRammler::sample() const
{
    scalar K = 1.0 - Foam::exp(-Foam::pow((maxValue_ - minValue_)/d_, n_));
    scalar y = rndGen_.sample01<scalar>();
    scalar x = minValue_ + d_*Foam::pow(-Foam::log(1.0 - y*K), 1.0/n_);
    return x;
}


Foam::scalar Foam::distributionModels::RosinRammler::minValue() const
{
    return minValue_;
}


Foam::scalar Foam::distributionModels::RosinRammler::maxValue() const
{
    return maxValue_;
}


Foam::scalar Foam::distributionModels::RosinRammler::meanValue() const
{
    return d_;
}


// ************************************************************************* //
