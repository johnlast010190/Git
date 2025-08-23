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
    (c) held by original author

\*---------------------------------------------------------------------------*/

#include "potentialCurrent.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(potentialCurrent, 0);
addToRunTimeSelectionTable(waveTheory, potentialCurrent, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


potentialCurrent::potentialCurrent
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    waveTheory(subDictName, mesh_),
    U_(vector(coeffDict_.lookup("U"))),
    Tsoft_(coeffDict_.lookup<scalar>("Tsoft")),
    localSeaLevel_
    (
        coeffDict_.lookupOrDefault<scalar>("localSeaLevel", seaLevel_)
    )
{}


void potentialCurrent::printCoeffs()
{
    Info<< "Loading wave theory: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar potentialCurrent::factor(const scalar& time) const
{
    scalar factor(1);
    if (Tsoft_ > 0.0)
    {
        factor = Foam::sin(PI_/2.0/Tsoft_*Foam::min(Tsoft_, time));
    }

    return factor;
}


scalar potentialCurrent::eta
(
    const point& x,
    const scalar& time
) const
{
//    scalar eta = seaLevel_;
    scalar eta = localSeaLevel_;
    return eta;
}


scalar potentialCurrent::ddxPd
(
    const point& x,
    const scalar& time,
    const vector& unitVector
) const
{
    return 0.0;
}


scalar potentialCurrent::p
(
    const point& x,
    const scalar& time
) const
{
    scalar result = rhoWater_*Foam::mag(g_)*localSeaLevel_;

    return result;
}


vector potentialCurrent::U
(
    const point& x,
    const scalar& time
) const
{
    return (U_*factor(time));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
