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
    (c) 2018 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "coordinate/rotation/cylindricalRotation.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace coordinateRotations
    {
        defineTypeName(cylindrical);
        addToRunTimeSelectionTable
        (
            coordinateRotation,
            cylindrical,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::tensor Foam::coordinateRotations::cylindrical::rotation
(
    const vector& axis
)
{
    // Guess second axis
    return axes::rotation(axis, Zero, E3_E1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateRotations::cylindrical::cylindrical(const cylindrical& crot)
:
    coordinateRotations::axes(crot)
{}


Foam::coordinateRotations::cylindrical::cylindrical(const vector& axis)
:
    coordinateRotations::axes(axis)  // Guess second axis
{}


Foam::coordinateRotations::cylindrical::cylindrical(const dictionary& dict)
:
    cylindrical(dict.getCompat<vector>("axis", {{"e3", -1806}}))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::coordinateRotations::cylindrical::write(Ostream& os) const
{
    os << type() << " axis: " << axis1_;
}


void Foam::coordinateRotations::cylindrical::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    os.beginBlock(keyword);

    os.writeEntry("type", type());
    os.writeEntry("axis", axis1_);

    os.endBlock();
}


// ************************************************************************* //
