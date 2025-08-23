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

Class
    Foam::hPolynomialThermo

Group
    grpSpecieThermo

Description
    Enthalpy based thermodynamics package using a polynomial function
    of temperature for the constant heat capacity at constant pressure:

    \verbatim
        Cp = a0 + a1*T + a2*sqr(T) + a3*pow3(T) + a4*pow4(T)
           + a5*pow(T, 5) + a6*pow(T, 6) + a7*pow(T, 7)
    \endverbatim

    The polynomial function is templated on the order of the polynomial which
    defaults to 8.

Usage

    \table
        Property     | Description
        Hf           | Heat of formation
        Sf           | Standard entropy
        CpCoeffs<8>  | Specific heat at constant pressure polynomial coeffs
    \endtable

    Example of the specification of the thermodynamic properties:
    \verbatim
    thermodynamics
    {
        Hf              0;
        Sf              0;
        CpCoeffs<8>     (1000 -0.05 0.003 0 0 0 0 0);
    }
    \endverbatim


Note
    Specific heat at constant pressure polynomial coefficients evaluate to an
    expression in [J/kg/K].


SourceFiles
    hPolynomialThermoI.H
    hPolynomialThermo.C

See also
    Foam::Polynomial

\*---------------------------------------------------------------------------*/

#ifndef hPolynomialThermo_H
#define hPolynomialThermo_H

#include "primitives/Scalar/scalar/scalar.H"
#include "primitives/functions/Polynomial/Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState, int PolySize>
class hPolynomialThermo;

template<class EquationOfState, int PolySize>
inline hPolynomialThermo<EquationOfState, PolySize> operator+
(
    const hPolynomialThermo<EquationOfState, PolySize>&,
    const hPolynomialThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
inline hPolynomialThermo<EquationOfState, PolySize> operator*
(
    const scalar,
    const hPolynomialThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
inline hPolynomialThermo<EquationOfState, PolySize> operator==
(
    const hPolynomialThermo<EquationOfState, PolySize>&,
    const hPolynomialThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const hPolynomialThermo<EquationOfState, PolySize>&
);

class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class hPolynomialThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState, int PolySize=8>
class hPolynomialThermo
:
    public EquationOfState
{
    // Private data

        //- Heat of formation
        scalar hf_;

        //- Standard entropy
        scalar Sf_;

        //- Specific heat at constant pressure polynomial coeffs [J/kg/K/K^i]
        Polynomial<PolySize> CpCoeffs_;

        //- Enthalpy polynomial coeffs [J/kg/K^i]
        //  Derived from Cp coeffs. Relative to Tstd.
        typename Polynomial<PolySize>::intPolyType hCoeffs_;

        //- Entropy polynomial coeffs [J/kg/K/K^i]
        //  Derived from Cp coeffs. Relative to Tstd.
        Polynomial<PolySize> sCoeffs_;


    // Private Member Functions

        //- Construct from components
        inline hPolynomialThermo
        (
            const EquationOfState& pt,
            const scalar hf,
            const scalar Sf,
            const Polynomial<PolySize>& CpCoeffs,
            const typename Polynomial<PolySize>::intPolyType& hCoeffs,
            const Polynomial<PolySize>& sCoeffs
        );


public:

    // Constructors

        //- Construct as copy
        hPolynomialThermo(const hPolynomialThermo&) = default;

        //- Construct from dictionary
        hPolynomialThermo(const objectRegistry& obr, const dictionary& dict);

        //- Construct as a named copy
        inline hPolynomialThermo(const word&, const hPolynomialThermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hPolynomial<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar) const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/kg/K]
            inline scalar Cp(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kg]
            inline scalar ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar hs(const scalar p, const scalar T) const;

            //- Enthalpy of formation [J/kg]
            inline scalar hf() const;

            //- Entropy [J/kg/K]
            inline scalar s(const scalar p, const scalar T) const;

            //- Gibbs free energy of the mixture in the standard state [J/kg]
            inline scalar gStd(const scalar T) const;

            //- Including conversion functions from H to E
            #include "../thermo/HtoEthermo.H"


        // Derivative term used for Jacobian

            //- Temperature derivative of heat capacity at constant pressure
            inline scalar dCpdT(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator=(const hPolynomialThermo&);
        inline void operator+=(const hPolynomialThermo&);
        inline void operator*=(const scalar);


    // Friend operators

        friend hPolynomialThermo operator+ <EquationOfState, PolySize>
        (
            const hPolynomialThermo&,
            const hPolynomialThermo&
        );

        friend hPolynomialThermo operator* <EquationOfState, PolySize>
        (
            const scalar,
            const hPolynomialThermo&
        );

        friend hPolynomialThermo operator== <EquationOfState, PolySize>
        (
            const hPolynomialThermo&,
            const hPolynomialThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState, PolySize>
        (
            Ostream&,
            const hPolynomialThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermo/hPolynomial/hPolynomialThermoI.H"

#ifdef NoRepository
    #include "thermo/hPolynomial/hPolynomialThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
