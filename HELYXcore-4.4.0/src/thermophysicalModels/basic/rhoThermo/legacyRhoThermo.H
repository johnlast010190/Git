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
    (c) 2015-2017 OpenCFD Ltd.

Class
    Foam::legacyRhoThermo

Description
    Energy for a mixture based on density

SourceFiles
    legacyRhoThermo.C

\*---------------------------------------------------------------------------*/

#ifndef legacyRhoThermo_H
#define legacyRhoThermo_H

#include "rhoThermo/rhoThermo.H"
#include "basicThermo/legacyBasicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class legacyRhoThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicRhoThermo, class MixtureType>
class legacyRhoThermo
:
    public legacyBasicThermo<BasicRhoThermo, MixtureType>
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate
        (
            const volScalarField& p,
            volScalarField& T,
            volScalarField& he,
            volScalarField& psi,
            volScalarField& rho,
            volScalarField& mu,
            volScalarField& alpha,
            volScalarField& Cp,
            volScalarField& Cv,
            const bool doOldTimes
        );

        //- Construct as copy (not implemented)
        legacyRhoThermo(const legacyRhoThermo<BasicRhoThermo, MixtureType>&);


public:

    //- Runtime type information
    TypeName("heRhoThermo");


    // Constructors

        //- Construct from obr and phase name
        legacyRhoThermo(const objectRegistry& obr, const word& phaseName);


    //- Destructor
    virtual ~legacyRhoThermo();


    // Member functions

        //- Update properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "rhoThermo/legacyRhoThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
