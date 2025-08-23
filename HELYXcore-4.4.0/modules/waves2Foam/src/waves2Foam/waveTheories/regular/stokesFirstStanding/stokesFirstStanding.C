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

#include "stokesFirstStanding.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(stokesFirstStanding, 0);
addToRunTimeSelectionTable(waveTheory, stokesFirstStanding, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


stokesFirstStanding::stokesFirstStanding
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    waveTheory(subDictName, mesh_),
    H_(coeffDict_.lookup<scalar>("height")),
    h_(coeffDict_.lookup<scalar>("depth")),
    omega_(coeffDict_.lookup<scalar>("omega")),
    period_(2*PI_/omega_),
    phi_(coeffDict_.lookup<scalar>("phi")),
    k_(vector(coeffDict_.lookup("waveNumber"))),
    K_(mag(k_)),

    Tsoft_(coeffDict_.lookupOrDefault<scalar>("Tsoft",period_))
{
    checkWaveDirection(k_);
}


void stokesFirstStanding::printCoeffs()
{
    Info<< "Loading wave theory: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar stokesFirstStanding::factor(const scalar& time) const
{
    scalar factor(1.0);
    if (Tsoft_ > 0.0)
    {
        factor = Foam::sin(omega_/4.0*Foam::min(Tsoft_, time));
    }

    return factor;
}


scalar stokesFirstStanding::eta
(
    const point& x,
    const scalar& time
) const
{
    scalar eta =
        (
            H_/2.0*Foam::cos(omega_*time - (k_ & x) + phi_)
          + H_/2.0*Foam::cos(omega_*time + (k_ & x) + phi_)
        )*factor(time) + seaLevel_;
    return eta;
}


scalar stokesFirstStanding::ddxPd
(
    const point& x,
    const scalar& time,
    const vector& unitVector
) const
{

    scalar Z(returnZ(x));
    scalar arg1(omega_*time - (k_ & x) + phi_);
    scalar arg2(omega_*time + (k_ & x) + phi_);

    scalar ddxPd(0);

    ddxPd = (
                rhoWater_*mag(g_)*K_*H_/2.0*Foam::cosh(K_*(Z + h_))
               /Foam::cosh(K_*h_)*Foam::sin(arg1)
              - rhoWater_*mag(g_)*K_*H_/2.0*Foam::cosh(K_*(Z + h_))
               /Foam::cosh(K_*h_)*Foam::sin(arg2)
            )*factor(time);

    return ddxPd;
}


vector stokesFirstStanding::U
(
    const point& x,
    const scalar& time
) const
{
    scalar Z(returnZ(x));

    scalar Uhorz = PI_*H_/period_ *
                   Foam::cosh(K_*(Z + h_))/Foam::sinh(K_*h_) *
                   Foam::cos(omega_*time - (k_ & x) + phi_)
                 - PI_*H_/period_ *
                   Foam::cosh(K_*(Z + h_))/Foam::sinh(K_*h_) *
                   Foam::cos(omega_*time + (k_ & x) + phi_);

    Uhorz *= factor(time);

    scalar Uvert = - PI_*H_/period_ *
                   Foam::sinh(K_*(Z + h_))/Foam::sinh(K_*h_) *
                   Foam::sin(omega_*time - (k_ & x) + phi_)
                 - PI_*H_/period_ *
                   Foam::sinh(K_*(Z + h_))/Foam::sinh(K_*h_) *
                   Foam::sin(omega_*time + (k_ & x) + phi_);

    Uvert *= factor(time);

    // Note "-" because of "g" working in the opposite direction
    return Uhorz*k_/K_ - Uvert*direction_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
