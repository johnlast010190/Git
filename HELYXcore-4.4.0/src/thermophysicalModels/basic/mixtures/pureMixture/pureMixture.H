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
    (c) 2011-2020 OpenFOAM Foundation

Class
    Foam::pureMixture

Description
    Foam::pureMixture

SourceFiles
    pureMixture.C

\*---------------------------------------------------------------------------*/

#ifndef pureMixture_H
#define pureMixture_H

#include "mixtures/basicMixture/basicMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pureMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class pureMixture
:
    public basicMixture
{
public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;

    //- Mixing type for thermodynamic properties
    typedef ThermoType thermoMixtureType;

    //- Mixing type for transport properties
    typedef ThermoType transportMixtureType;


private:

    // Private Data

        ThermoType mixture_;


public:
    // Constructors

        //- Construct from dictionary, mesh and phase name
        pureMixture(const dictionary&, const objectRegistry&, const word&);

        //- Disallow default bitwise copy construction
        pureMixture(const pureMixture<ThermoType>&) = delete;


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "pureMixture<" + ThermoType::typeName() + '>';
        }

        const thermoMixtureType& cellThermoMixture(const label) const
        {
            return mixture_;
        }

        const thermoMixtureType& patchFaceThermoMixture
        (
            const label,
            const label
        ) const
        {
            return mixture_;
        }

        const transportMixtureType& cellTransportMixture(const label) const
        {
            return mixture_;
        }

        const transportMixtureType& patchFaceTransportMixture
        (
            const label,
            const label
        ) const
        {
            return mixture_;
        }

        const transportMixtureType& cellTransportMixture
        (
            const label,
            const thermoMixtureType&
        ) const
        {
            return mixture_;
        }

        const transportMixtureType& patchFaceTransportMixture
        (
            const label,
            const label,
            const thermoMixtureType&
        ) const
        {
            return mixture_;
        }

        //- Read dictionary
        void read(const objectRegistry&, const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "mixtures/pureMixture/pureMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
