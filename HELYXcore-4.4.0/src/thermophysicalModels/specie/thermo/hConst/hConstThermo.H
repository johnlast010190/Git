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
    Foam::hConstThermo

Group
    grpSpecieThermo

Description
    Enthalpy based thermodynamics package using a constant heat capacity
    at constant pressure:

    \verbatim
        h = Cp*(T - Tref) + Hsref
    \endverbatim

Usage
    \table
        Property     | Description
        Cp           | Constant Heat capacity at constant pressure [J/kg/K]
        Tref         | Reference temperature [K] (defaults to Tstd)
        Hsref        | Reference sensible enthalpy [J/kg] (defaults to 0)
        Hf           | Heat of formation [J/kg]
    \endtable

    Example specification of hConstThermo for air:

    \verbatim
    thermodynamics
    {
        Cp          1007;
        Hf          0;
    }
    \endverbatim

SourceFiles
    hConstThermoI.H
    hConstThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hConstThermo_H
#define hConstThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class hConstThermo;

template<class EquationOfState>
inline hConstThermo<EquationOfState> operator+
(
    const hConstThermo<EquationOfState>&,
    const hConstThermo<EquationOfState>&
);

template<class EquationOfState>
inline hConstThermo<EquationOfState> operator*
(
    const scalar,
    const hConstThermo<EquationOfState>&
);

template<class EquationOfState>
inline hConstThermo<EquationOfState> operator==
(
    const hConstThermo<EquationOfState>&,
    const hConstThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const hConstThermo<EquationOfState>&
);

class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class hConstThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hConstThermo
:
    public EquationOfState
{
    // Private data

        //- Heat capacity at constant pressure [J/kg/K]
        scalar Cp_;

        //- Heat of formation [J/kg]
        scalar hf_;

        //- Reference temperature around which to linearise [K]
        scalar Tref_;

        //- Reference sensible enthalpy around which to linearise [J/kg]
        scalar hsref_;

        scalar Sf_;


    // Private Member Functions

        //- Construct from components
        inline hConstThermo
        (
            const EquationOfState& st,
            const scalar Cp,
            const scalar hf,
            const scalar Tref,
            const scalar hsref,
            const scalar Sf
        );


public:

    // Constructors

        //- Construct from dictionary
        hConstThermo(const objectRegistry& obr, const dictionary& dict);

        //- Construct as named copy
        inline hConstThermo(const word&, const hConstThermo&);

        //- Construct and return a clone
        inline autoPtr<hConstThermo> clone() const;

        //- Selector from dictionary
        inline static autoPtr<hConstThermo> New
        (
            const objectRegistry& obr,
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hConst<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/kg/K]
            inline scalar Cp(const scalar p, const scalar T) const;

            //- Absolute enthalpy [J/kg]
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

        inline void operator+=(const hConstThermo&);


    // Friend operators

        friend hConstThermo operator+ <EquationOfState>
        (
            const hConstThermo&,
            const hConstThermo&
        );

        friend hConstThermo operator* <EquationOfState>
        (
            const scalar,
            const hConstThermo&
        );

        friend hConstThermo operator== <EquationOfState>
        (
            const hConstThermo&,
            const hConstThermo&
        );


    // IOstream Operators

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const hConstThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermo/hConst/hConstThermoI.H"

#ifdef NoRepository
    #include "thermo/hConst/hConstThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
