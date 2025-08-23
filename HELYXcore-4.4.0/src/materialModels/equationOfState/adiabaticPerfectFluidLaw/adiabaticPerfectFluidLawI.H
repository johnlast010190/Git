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
    (c) 2013-2017 OpenFOAM Foundation
    (c) 2016-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "adiabaticPerfectFluidLaw.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::adiabaticPerfectFluidLaw>
Foam::adiabaticPerfectFluidLaw::clone() const
{
    return autoPtr<adiabaticPerfectFluidLaw>
    (
        new adiabaticPerfectFluidLaw(*this)
    );
}


inline Foam::autoPtr<Foam::adiabaticPerfectFluidLaw>
Foam::adiabaticPerfectFluidLaw::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<adiabaticPerfectFluidLaw>
    (
        new adiabaticPerfectFluidLaw(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::adiabaticPerfectFluidLaw::rho(const scalar p) const
{
    return rho0_*pow((p + B_)/(p0_ + B_), 1.0/gamma_);
}


inline Foam::scalar Foam::adiabaticPerfectFluidLaw::spContribution
(
    const scalar p,
    const scalar T
) const
{
    const scalar n = 1 - 1.0/gamma_;
    return
       -pow(p0_ + B_, 1.0/gamma_)
       *(pow((p + B_), n) - pow((Pstd + B_), n))
       /(rho0_*T*n);
}


inline Foam::scalar Foam::adiabaticPerfectFluidLaw::svContribution() const
{
    NotImplemented;
}


inline Foam::scalar Foam::adiabaticPerfectFluidLaw::Z
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    return p/(rho(p)*R*T);
}


inline Foam::scalar Foam::adiabaticPerfectFluidLaw::psi(const scalar p) const
{
    return
        (rho0_/(gamma_*(p0_ + B_)))
       *pow((p + B_)/(p0_ + B_), 1.0/gamma_ - 1.0);
}


// ************************************************************************* //
