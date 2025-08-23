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

#include "templateWaveTheory.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(templateWaveTheory, 0);
addToRunTimeSelectionTable(waveTheory, templateWaveTheory, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

templateWaveTheory::templateWaveTheory
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    waveTheory(subDictName, mesh_)//,
    //H_(coeffDict_.lookup<scalar>("height")),
    //Add the other variables
    //Tsoft_(coeffDict_.lookupOrDefault<scalar>("Tsoft",period_))
{}


void templateWaveTheory::printCoeffs()
{
    Info<< "Loading wave theory: " << typeName << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar templateWaveTheory::factor(const scalar& time) const
{
    scalar factor(1.0);
    if (Tsoft_ > 0.0)
    {
        factor = Foam::sin(2*PI_/(4.0*Tsoft_)*Foam::min(Tsoft_, time));
    }

    return factor;
}

scalar templateWaveTheory::eta
(
    const point& x,
    const scalar& time
) const
{
    scalar eta = 0.0;

    // Insert expression for eta

    return eta;
}

scalar templateWaveTheory::ddxPd
(
    const point& x,
    const scalar& time,
    const vector& unitVector
) const
{
    // This gives the z-coordinate relative to seaLevel
    scalar Z(returnZ(x));

    scalar ddxPd(0);
    // Most theories return 0, as the issue with oblique waves has not
    // yet been derived.

    return ddxPd;
}

vector templateWaveTheory::U
(
    const point& x,
    const scalar& time
) const
{
    // This gives the z-coordinate relative to seaLevel
    scalar Z(returnZ(x));

    scalar Uhorz(0.0), Uvert(0);

    // Insert expression for Uvert and Uhorz as a function of z-coordinate

    // Scale by ramp-up time
    Uhorz *= factor(time);
    Uvert *= factor(time);

    // Cast into a directed vector. Look into any of the existing waves
    // theories to test how this is done correctly. Without the definitions
    // of the variables, it looks something like:
    //  return Uhorz*k_/K_ - Uvert*direction_;
    // Note "-" because of "g" working in the opposite direction

    return vector::zero;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
