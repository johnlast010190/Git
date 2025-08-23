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
    Foam::sensibleInternalEnergy

Group
    grpSpecieThermo

Description
    Thermodynamics mapping class to expose the sensible internal energy
    functions.

\*---------------------------------------------------------------------------*/

#ifndef sensibleInternalEnergy_H
#define sensibleInternalEnergy_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sensibleInternalEnergy Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class sensibleInternalEnergy
{

public:

    // Constructors

        //- Construct
        sensibleInternalEnergy()
        {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "sensibleInternalEnergy";
        }

        // Fundamental properties

            static bool enthalpy()
            {
                return false;
            }

            static word energyName()
            {
                return "e";
            }

            //- Heat capacity at constant volume [J/kg/K]
            scalar Cpv
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                return thermo.Cv(p, T);
            }

            //- Sensible internal energy [J/kg]
            scalar he
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                return thermo.es(p, T);
            }

            //- Temperature from sensible internal energy
            //  given an initial temperature T0
            scalar The
            (
                const Thermo& thermo,
                const scalar e,
                const scalar p,
                const scalar T0
            ) const
            {
                return thermo.Tes(e, p, T0);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
