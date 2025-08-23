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
    (c) 2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "massRosinRammler/massRosinRammler.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{
    defineTypeNameAndDebug(massRosinRammler, 0);
    addToRunTimeSelectionTable(distributionModel, massRosinRammler, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionModels::massRosinRammler::massRosinRammler
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


Foam::distributionModels::massRosinRammler::massRosinRammler
(
    const massRosinRammler& p
)
:
    distributionModel(p),
    minValue_(p.minValue_),
    maxValue_(p.maxValue_),
    d_(p.d_),
    n_(p.n_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributionModels::massRosinRammler::~massRosinRammler()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::massRosinRammler::sample() const
{
    scalar d;

    // Re-sample if the calculated d is out of the physical range
    do
    {
        const scalar a = 3/n_ + 1;
        const scalar P = rndGen_.sample01<scalar>();
        const scalar x = invIncGamma(a, P);
        d = d_*pow(x, 1/n_);
    } while (d < minValue_ || d > maxValue_);

    return d;
}


Foam::scalar Foam::distributionModels::massRosinRammler::minValue() const
{
    return minValue_;
}


Foam::scalar Foam::distributionModels::massRosinRammler::maxValue() const
{
    return maxValue_;
}


Foam::scalar Foam::distributionModels::massRosinRammler::meanValue() const
{
    return d_;
}


// ************************************************************************* //
