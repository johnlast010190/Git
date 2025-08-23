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

#include "linear.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace partitioningModels
{
    defineTypeNameAndDebug(linear, 0);
    addToRunTimeSelectionTable
    (
        partitioningModel,
        linear,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::partitioningModels::
linear::linear(const dictionary& dict)
:
    partitioningModel(),
    alphaLiquid1_(dict.lookup<scalar>("alphaLiquid1")),
    alphaLiquid0_(dict.lookup<scalar>("alphaLiquid0"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoilingModels::partitioningModels::
linear::~linear()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::partitioningModels::
linear::fLiquid
(
    const scalarField& alphaLiquid
) const
{
    return max
        (
            scalar(0),
            min
            (
                scalar(1)-(alphaLiquid1_ - alphaLiquid)
               /(alphaLiquid1_ - alphaLiquid0_),
                scalar(1)
            )
        );
}


void Foam::wallBoilingModels::partitioningModels::
linear::write(Ostream& os) const
{
    partitioningModel::write(os);
    os.writeEntry("alphaLiquid1", alphaLiquid1_);
    os.writeEntry("alphaLiquid0", alphaLiquid0_);
}


// ************************************************************************* //
