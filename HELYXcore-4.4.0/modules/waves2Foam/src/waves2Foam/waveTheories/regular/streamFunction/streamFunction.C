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

#include "streamFunction.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(streamFunction, 0);
addToRunTimeSelectionTable(waveTheory, streamFunction, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


streamFunction::streamFunction
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    waveTheory(subDictName, mesh_),
    h_(coeffDict_.lookup<scalar>("depth")),
    omega_(coeffDict_.lookup<scalar>("omega")),
    period_(2*PI_/omega_),
    phi_(coeffDict_.lookup<scalar>("phi")),
    k_(vector(coeffDict_.lookup("waveNumber"))),
    K_(mag(k_)),
    uBar_(coeffDict_.lookup<scalar>("uBar")),
    N_(coeffDict_.lookup<scalar>("N")),
    A_("A", coeffDict_, N_),
    B_("B", coeffDict_, N_),

    Tsoft_(coeffDict_.lookupOrDefault<scalar>("Tsoft",period_)),
    Tstart_(coeffDict_.lookupOrDefault<scalar>("Tstart",0.0)),
    Tend_(coeffDict_.lookupOrDefault<scalar>("Tend",GREAT))
{
    checkWaveDirection(k_);
}


void streamFunction::printCoeffs()
{
    Info<< "Loading wave theory: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar streamFunction::factor(const scalar& time) const
{
    scalar factor(1.0);
    if (time < Tstart_)
    {
        factor = 0.0;
    }
    else if (Tsoft_ > 0.0)
    {
        factor = Foam::sin(2*PI_/(4.0*Tsoft_)
            *Foam::min(Tsoft_, time - Tstart_));
    }

    if (time > Tend_)
    {
        factor = Foam::cos(2*PI_/(4.0*Tsoft_)*Foam::min(Tsoft_, time - Tend_));
    }

    return factor;
}


scalar streamFunction::eta
(
    const point& x,
    const scalar& time
) const
{
    scalar eta(0);
    scalar arg = (k_ & x) - omega_*time + phi_;

    forAll(A_,ii)
    {
        eta += A_[ii]*Foam::cos((ii + 1)*arg);
    }
    eta *= factor(time);
    eta += seaLevel_;
    return eta;
}


scalar streamFunction::ddxPd
(
    const point& x,
    const scalar& time,
    const vector& unitVector
) const
{
    // PAS PAA, DETTE FOELGER IKKE KONVENTIONEN MED ARBITRAER BOELGETALSVEKTOR!
    scalar Z(returnZ(x));
    scalar arg = (k_ & x) - omega_*time + phi_;

    scalar ddxPd(0);
    vector uu(U(x,time));
    scalar Ux(uu.x());
    scalar Uy(- uu & direction_);
    scalar dUx(0);
    scalar dUy(0);

    scalar j(0);
    forAll(B_,ii)
    {
        j = ii + 1;
        dUx += - j*K_*B_[ii]*Foam::cosh(j*K_*(Z + h_))/Foam::cosh(j*K_*h_)
            *Foam::sin(j*arg);
        dUy +=   j*K_*B_[ii]*Foam::sinh(j*K_*(Z + h_))/Foam::cosh(j*K_*h_)
            *Foam::cos(j*arg);
    }
    dUx *= factor(time);
    dUy *= factor(time);

    ddxPd =  rhoWater_*(omega_/K_*dUx - Ux*dUx - Uy*dUy)*factor(time);

    return ddxPd;
}


vector streamFunction::U
(
    const point& x,
    const scalar& time
) const
{
    scalar Z(returnZ(x));
    scalar Uhorz(omega_/K_ - uBar_);
    scalar Uvert(0);
    scalar arg = (k_ & x) - omega_*time + phi_;
    scalar j(0);

    forAll(B_,ii)
    {
        j = ii + 1;

        Uhorz += B_[ii]*Foam::cosh(j*K_*(Z + h_))/Foam::cosh(j*K_*h_)
            *Foam::cos(j*arg);

        Uvert += B_[ii]*Foam::sinh(j*K_*(Z + h_))/Foam::cosh(j*K_*h_)
            *Foam::sin(j*arg);
    }
    Uhorz *= factor(time);
    Uvert *= factor(time);

    // Note "-" because of "g" working in the opposite direction
    return Uhorz*k_/K_ - Uvert*direction_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
