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
    (c) 2011-2023 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

inline Foam::scalar Foam::liquidProperties::limit(const scalar T) const
{
    return T;
}


inline Foam::scalar Foam::liquidProperties::Y() const
{
    return 1;
}


inline Foam::scalar Foam::liquidProperties::Tc() const
{
    return Tc_;
}


inline Foam::scalar Foam::liquidProperties::Pc() const
{
    return Pc_;
}


inline Foam::scalar Foam::liquidProperties::Vc() const
{
    return Vc_;
}


inline Foam::scalar Foam::liquidProperties::Zc() const
{
    return Zc_;
}


inline Foam::scalar Foam::liquidProperties::Tt() const
{
    return Tt_;
}


inline Foam::scalar Foam::liquidProperties::Pt() const
{
    return Pt_;
}


inline Foam::scalar Foam::liquidProperties::Tb() const
{
    return Tb_;
}


inline Foam::scalar Foam::liquidProperties::dipm() const
{
    return dipm_;
}


inline Foam::scalar Foam::liquidProperties::omega() const
{
    return omega_;
}


inline Foam::scalar Foam::liquidProperties::delta() const
{
    return delta_;
}


inline Foam::scalar Foam::liquidProperties::psi(scalar p, scalar T) const
{
    return 0;
}


inline Foam::scalar Foam::liquidProperties::CpMCv(scalar p, scalar T) const
{
    return 0;
}


inline Foam::scalar Foam::liquidProperties::es(scalar p, const scalar T) const
{
    return hs(p, T) - p/rho(p, T);
}


inline Foam::scalar Foam::liquidProperties::ea(scalar p, const scalar T) const
{
    return ha(p, T) - p/rho(p, T);
}


inline Foam::scalar Foam::liquidProperties::rho(scalar p, scalar T) const
{
    return rho_().value(T);
}


inline Foam::scalar Foam::liquidProperties::pv(scalar p, scalar T) const
{
    return pv_().value(T);
}


inline Foam::scalar Foam::liquidProperties::hl(scalar p, scalar T) const
{
    return hl_().value(T);
}


inline Foam::scalar Foam::liquidProperties::Cp(scalar p, scalar T) const
{
    return CpContribution(p, T) + Cp_().value(T);
}


inline Foam::scalar Foam::liquidProperties::hContribution
(
    scalar p,
    scalar T
) const
{
    return (p - pRef_)/rho(p, T);
}


inline Foam::scalar Foam::liquidProperties::CpContribution
(
    scalar p,
    scalar T
) const
{
    return -((p - pRef_)/sqr(rho(p, T)))*rho_->derivative(T);
}


inline Foam::scalar Foam::liquidProperties::hs(scalar p, scalar T) const
{
    return ha(p, T) - hf();
}


inline Foam::scalar Foam::liquidProperties::hf() const
{
    return hf_;
}


inline Foam::scalar Foam::liquidProperties::ha(scalar p, scalar T) const
{
    return hContribution(p, T) + h_().value(T);
}


inline Foam::scalar Foam::liquidProperties::Cpg(scalar p, scalar T) const
{
    return Cpg_().value(T);
}


inline Foam::scalar Foam::liquidProperties::B(scalar p, scalar T) const
{
    return B_().value(T);
}


inline Foam::scalar Foam::liquidProperties::mu(scalar p, scalar T) const
{
    return mu_().value(T);
}


inline Foam::scalar Foam::liquidProperties::mug(scalar p, scalar T) const
{
    return mug_().value(T);
}


inline Foam::scalar Foam::liquidProperties::kappa(scalar p, scalar T) const
{
    return kappa_().value(T);
}


inline Foam::scalar Foam::liquidProperties::kappag(scalar p, scalar T) const
{
    return kappag_().value(T);
}


inline Foam::scalar Foam::liquidProperties::sigma(scalar p, scalar T) const
{
    return sigma_().value(T);
}


inline Foam::scalar Foam::liquidProperties::D(scalar p, scalar T) const
{
    return Dfunc1_.valid() ? Dfunc1_().value(T) : D_().value(p, T);
}


inline Foam::scalar Foam::liquidProperties::D(scalar p, scalar T, scalar Wb) const
{
    return Dfunc1_.valid() ? Dfunc1_().value(T) : D_().value(p, T, Wb);
}


// ************************************************************************* //
