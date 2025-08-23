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
    (c) 2012-2023 OpenFOAM Foundation

Class
    Foam::sensibleEnthalpy

Group
    grpSpecieThermo

Description
    Thermodynamics mapping class to expose the sensible enthalpy functions.

\*---------------------------------------------------------------------------*/

#ifndef sensibleEnthalpy_H
#define sensibleEnthalpy_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sensibleEnthalpy Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class sensibleEnthalpy
{

public:

    // Constructors

        //- Construct
        sensibleEnthalpy()
        {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "sensibleEnthalpy";
        }

        // Fundamental properties

            static bool enthalpy()
            {
                return true;
            }

            static word energyName()
            {
                return "h";
            }

            // Heat capacity at constant pressure [J/kg/K]
            scalar Cpv
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                return thermo.Cp(p, T);
            }

            // Sensible enthalpy [J/kg]
            scalar he
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                return thermo.hs(p, T);
            }

            //- Temperature from sensible enthalpy
            //  given an initial temperature T0
            scalar The
            (
                const Thermo& thermo,
                const scalar h,
                const scalar p,
                const scalar T0
            ) const
            {
                return thermo.Ths(h, p, T0);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
