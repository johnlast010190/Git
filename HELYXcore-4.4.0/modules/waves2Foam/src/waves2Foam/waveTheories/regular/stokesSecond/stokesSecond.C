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

#include "stokesSecond.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(stokesSecond, 0);
addToRunTimeSelectionTable(waveTheory, stokesSecond, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


stokesSecond::stokesSecond
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
    Tsoft_(coeffDict_.lookupOrDefault<scalar>("Tsoft",period_)),
    debug_(Switch(coeffDict_.lookup("debug")))
{
    checkWaveDirection(k_);

    if
    (
        H_/2.0 - 4.0*1.0/16.0*K_*sqr(H_)*(3.0/Foam::pow(Foam::tanh(K_*h_),3.0)
        - 1.0/Foam::tanh(K_*h_)) < 0
    )
    {
        if (debug_)
        {
            WarningIn
            (
                "label stokesSecond::eta(point x, scalar time)"
            ) << endl << "The validity of stokes second order is violated."
            << endl << "a_1 < 4 a_2, being first and second order"
            << " amplitudes respectively." << endl << endl;
            Info<< "a1 = " << H_/2.0 << " , a2 = "
                 << (1.0/16.0*K_*sqr(H_)*(3.0/Foam::pow(Foam::tanh(K_*h_),3.0)
                    - 1.0/Foam::tanh(K_*h_))) << endl;
        }
    }
}


void stokesSecond::printCoeffs()
{
    Info<< "Loading wave theory: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar stokesSecond::factor(const scalar& time) const
{
    scalar factor(1.0);
    if (Tsoft_ > 0.0)
    {
        factor = Foam::sin(2*PI_/(4.0*Tsoft_)*Foam::min(Tsoft_, time));
    }

    return factor;
}


scalar stokesSecond::eta
(
    const point& x,
    const scalar& time
) const
{
    scalar arg(omega_*time - (k_ & x) + phi_);

    scalar eta = (
                     H_/2.0*Foam::cos(arg) // First order contribution.
                   // Second order contribution.
                   + 1.0/16.0*K_*sqr(H_)
                   *(
                        3.0/Foam::pow(Foam::tanh(K_*h_),3.0)
                      - 1.0/Foam::tanh(K_*h_)
                    )
                   *Foam::cos(2.0*arg)
                 )*factor(time)  // Hot-starting.
                 + seaLevel_;      // Adding sea level.

    return eta;
}


scalar stokesSecond::ddxPd
(
    const point& x,
    const scalar& time,
    const vector& unitVector
) const
{
    scalar Z(returnZ(x));
    scalar arg(omega_*time - (k_ & x) + phi_);

    scalar ddxPd(0);

    ddxPd = (
                rhoWater_*mag(g_)*K_*H_/2.0*Foam::cosh(K_*(Z + h_))
                /Foam::cosh(K_*h_)*Foam::sin(arg)
                + (scalar(1)/4)*rhoWater_*mag(g_)*pow(K_,2)*pow(H_,2)
                /Foam::sinh(2*K_*h_)
                *( 3*Foam::cosh(2*K_*(Z + h_))/pow(Foam::sinh(K_*h_),2) - 1)
                *Foam::sin(2*arg)
            )*factor(time);

    return ddxPd;
}


vector stokesSecond::U
(
    const point& x,
    const scalar& time
) const
{
    scalar Z(returnZ(x));
    scalar cel(omega_/K_);
    scalar arg(omega_*time - (k_ & x) + phi_);

    // First order contribution
    scalar Uhorz = PI_*H_/period_ *
                   Foam::cosh(K_*(Z + h_))/Foam::sinh(K_*h_) *
                   Foam::cos(arg);

    // Second order contribution
    Uhorz += 3.0/16.0*cel*Foam::sqr(K_*H_)*Foam::cosh(2*K_*(Z + h_))
            /Foam::pow(Foam::sinh(K_*h_),4.0)*Foam::cos(2*arg)
             - 1.0/8.0*mag(g_)*sqr(H_)/(cel*h_);

    // First order contribution
    scalar Uvert = - PI_*H_/period_ *
                   Foam::sinh(K_*(Z + h_))/Foam::sinh(K_*h_) *
                   Foam::sin(arg);

    // Second order contribution
    Uvert += - 3.0/16.0*cel*sqr(K_*H_)*Foam::sinh(2*K_*(Z + h_))
            /Foam::pow(Foam::sinh(K_*h_), 4.0)*Foam::sin(2*arg);

    // Multiply by the time stepping factor
    Uvert *= factor(time);
    Uhorz *= factor(time);

    // Generate the velocity vector
    // Note "-" because of "g" working in the opposite direction
    return Uhorz*k_/K_ - Uvert*direction_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
