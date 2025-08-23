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
    (c) 2020-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "thermo/ePolynomial/ePolynomialThermo.H"
#include "db/IOstreams/IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
Foam::ePolynomialThermo<EquationOfState, PolySize>::ePolynomialThermo
(
    const objectRegistry& obr,
    const dictionary& dict
)
:
    EquationOfState(obr, dict),
    hf_(dict.subDict("thermodynamics").lookup<scalar>("Hf")),
    Sf_(dict.subDict("thermodynamics").lookup<scalar>("Sf")),
    CvCoeffs_
    (
        dict.subDict("thermodynamics").lookup
        (
            "CvCoeffs<" + Foam::name(PolySize) + '>'
        )
    ),
    eCoeffs_(),
    sCoeffs_()
{
    eCoeffs_ = CvCoeffs_.integral();
    sCoeffs_ = CvCoeffs_.integralMinus1();

    // Offset e poly so that it is relative to the enthalpy at Tstd
    eCoeffs_[0] -= eCoeffs_.value(Tstd);

    // Offset s poly so that it is relative to the entropy at Tstd
    sCoeffs_[0] -= sCoeffs_.value(Tstd);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
void Foam::ePolynomialThermo<EquationOfState, PolySize>::write
(
    Ostream& os
) const
{
    EquationOfState::write(os);

    dictionary dict("thermodynamics");
    dict.add("Hf", hf_);
    dict.add("Sf", Sf_);
    dict.add
    (
        word("CvCoeffs<" + Foam::name(PolySize) + '>'),
        CvCoeffs_
    );
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ePolynomialThermo<EquationOfState, PolySize>& pt
)
{
    pt.write(os);
    return os;
}


// ************************************************************************* //