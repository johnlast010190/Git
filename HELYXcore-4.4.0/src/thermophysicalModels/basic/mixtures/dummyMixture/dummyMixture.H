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
    (c) 2011-2022 OpenFOAM Foundation

Class
    Foam::dummyMixture

Description
    Foam::dummyMixture

SourceFiles
    dummyMixture.C

\*---------------------------------------------------------------------------*/

#ifndef dummyMixture_H
#define dummyMixture_H

#include "mixtures/basicCombustionMixture/basicCombustionMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dummyMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class dummyMixture
:
    public basicCombustionMixture
{
public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;

    //- Mixing type for thermodynamic properties
    typedef ThermoType thermoMixtureType;

    //- Mixing type for transport properties
    typedef ThermoType transportMixtureType;


private:

    // Private data

        dimensionedScalar stoicRatio_;


public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        dummyMixture
        (
            const dictionary&,
            const objectRegistry&,
            const word&
        );

        //- Disallow default bitwise copy construction
        dummyMixture(const dummyMixture<ThermoType>&) = delete;


    //- Destructor
    virtual ~dummyMixture()
    {}


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "dummyMixture<" + ThermoType::typeName() + '>';
        }

        const dimensionedScalar& stoicRatio() const
        {
            return stoicRatio_;
        }

        //- Read dictionary
        void read(const objectRegistry&, const dictionary&);

        //- Return thermo based on index
        const ThermoType& specieThermo(const label speciei) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "mixtures/dummyMixture/dummyMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
