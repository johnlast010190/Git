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
    (c) 2014-2023 OpenFOAM Foundation

Class
    Foam::PengRobinsonGas

Group
    grpSpecieEquationOfState

Description
    PengRobinsonGas cubic equation of state for gases
    Coefficient mixing is very inaccurate and not supported, so this equation of
    state can only be applied to mixtures by using a thermophysical properties
    model which mixes specie property value.

Usage
    \table
        Property     | Description
        Tc           | Critical temperature [K]
        Vc           | Critical volume [m^3/kmol]
        Pc           | Critical pressure [Pa]
        omega        | Acentric factor [-]
    \endtable
    Example specification of the PengRobinsonGas equation of state for methane:
    \verbatim
    equationOfState
    {
        Tc        190.55;
        Vc        0.0285;
        Pc        4.595e6;
        omega     0.0115;
    }
    \endverbatim

SourceFiles
    PengRobinsonGasI.H
    PengRobinsonGas.C

\*---------------------------------------------------------------------------*/

#ifndef PengRobinsonGas_H
#define PengRobinsonGas_H

#include "memory/autoPtr/autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class PengRobinsonGas;

template<class Specie>
inline PengRobinsonGas<Specie> operator+
(
    const PengRobinsonGas<Specie>&,
    const PengRobinsonGas<Specie>&
);

template<class Specie>
inline PengRobinsonGas<Specie> operator*
(
    const scalar,
    const PengRobinsonGas<Specie>&
);

template<class Specie>
inline PengRobinsonGas<Specie> operator==
(
    const PengRobinsonGas<Specie>&,
    const PengRobinsonGas<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const PengRobinsonGas<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class PengRobinsonGas Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class PengRobinsonGas
:
    public Specie
{
    // Private data

        //- Critical Temperature [K]
        scalar Tc_;

        //- Critical volume [m^3/kmol]
        scalar Vc_;

        //- Critical compression factor [-]
        scalar Zc_;

        //- Critical Pressure [Pa]
        scalar Pc_;

        //- Acentric factor [-]
        scalar omega_;


public:

    // Constructors

        //- Construct from components
        inline PengRobinsonGas
        (
            const Specie& sp,
            const fvMesh& mesh,
            const scalar& Tc,
            const scalar& Vc,
            const scalar& Zc,
            const scalar& Pc,
            const scalar& omega
        );

        //- Construct from dictionary
        PengRobinsonGas(const objectRegistry& obr, const dictionary& dict);

        //- Construct as named copy
        inline PengRobinsonGas(const word& name, const PengRobinsonGas&);

        //- Construct and return a clone
        inline autoPtr<PengRobinsonGas> clone() const;

        // Selector from dictionary
        inline static autoPtr<PengRobinsonGas> New
        (
            const objectRegistry& obr,
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "PengRobinsonGas<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy contribution [J/kg]
            inline scalar h(const scalar p, const scalar T) const;

            //- Return Cp contribution [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return internal energy contribution [J/kg]
            inline scalar e(const scalar p, const scalar T) const;

            //- Return Cv contribution [J/(kg K]
            inline scalar Cv(scalar p, scalar T) const;

            //- Return entropy correction to the integral of Cp/T [J/kg/K]
            inline scalar sp(const scalar p, const scalar T) const;

            //- Return entropy correction to the integral of Cv/T [J/kg/K]
            inline scalar sv(const scalar p, const scalar T) const;

            //- Return compressibility [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;

            //- Return volumetric coefficient of thermal expansion [1/T]
            inline scalar alphav(const scalar p, const scalar T) const;


    // IO

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        //inline void operator=(const PengRobinsonGas&);
        inline void operator+=(const PengRobinsonGas&);
        inline void operator*=(const scalar);


    // Friend operators

        friend PengRobinsonGas operator+ <Specie>
        (
            const PengRobinsonGas&,
            const PengRobinsonGas&
        );

        friend PengRobinsonGas operator* <Specie>
        (
            const scalar s,
            const PengRobinsonGas&
        );

        friend PengRobinsonGas operator== <Specie>
        (
            const PengRobinsonGas&,
            const PengRobinsonGas&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const PengRobinsonGas&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "equationOfState/PengRobinsonGas/PengRobinsonGasI.H"

#ifdef NoRepository
    #include "equationOfState/PengRobinsonGas/PengRobinsonGas.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
