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
    Foam::absoluteInternalEnergy

Group
    grpSpecieThermo

Description
    Thermodynamics mapping class to expose the absolute internal energy
    functions.

\*---------------------------------------------------------------------------*/

#ifndef absoluteInternalEnergy_H
#define absoluteInternalEnergy_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class absoluteInternalEnergy Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class absoluteInternalEnergy
{

public:

    // Constructors

        //- Construct
        absoluteInternalEnergy()
        {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "absoluteInternalEnergy";
        }

        // Fundamental properties

            static bool enthalpy()
            {
                return false;
            }

            static word energyName()
            {
                return "ea";
            }

            // Heat capacity at constant volume [J/kg/K]
            scalar Cpv
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                return thermo.Cv(p, T);
            }

            // Absolute internal energy [J/kg]
            scalar he
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                return thermo.ea(p, T);
            }

            //- Temperature from absolute internal energy
            //  given an initial temperature T0
            scalar The
            (
                const Thermo& thermo,
                const scalar e,
                const scalar p,
                const scalar T0
            ) const
            {
                return thermo.Tea(e, p, T0);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
