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
    (c) 2012-2017 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "incompressiblePerfectGasLaw.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::incompressiblePerfectGasLaw>
Foam::incompressiblePerfectGasLaw::clone() const
{
    return autoPtr<incompressiblePerfectGasLaw>
    (
        new incompressiblePerfectGasLaw(*this)
    );
}


inline Foam::autoPtr<Foam::incompressiblePerfectGasLaw>
Foam::incompressiblePerfectGasLaw::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<incompressiblePerfectGasLaw>
    (
        new incompressiblePerfectGasLaw(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::incompressiblePerfectGasLaw::rho
(
    const scalar T,
    const scalar R
) const
{
    return pRef_/(R*T);
}


inline Foam::scalar Foam::incompressiblePerfectGasLaw::hContribution
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    return (p - Pstd)/rho(T, R);
}


inline Foam::scalar Foam::incompressiblePerfectGasLaw::CpContribution
(
    const scalar p,
    const scalar R
) const
{
    return (p - Pstd)/(pRef_/R);
}


inline Foam::scalar Foam::incompressiblePerfectGasLaw::Z
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    return p/(rho(T, R)*R*T);
}


inline Foam::scalar Foam::incompressiblePerfectGasLaw::alphav
(
    const scalar T
) const
{
    return 1/T;
}

// ************************************************************************* //
