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
    (c) 2015-2017 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::fv::buoyancyEnergy

Group
    grpFvOptionsSources

Description
    Calculates and applies the buoyancy energy source rho*(U&g) to the energy
    equation.

Usage
    Example usage:
    \verbatim
    fields          (h);                    // Name of energy field
    \endverbatim

SourceFiles
    buoyancyEnergy.C

\*---------------------------------------------------------------------------*/

#ifndef buoyancyEnergy_H
#define buoyancyEnergy_H

#include "cfdTools/general/fvOptions/fvOption.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class buoyancyEnergy Declaration
\*---------------------------------------------------------------------------*/

class buoyancyEnergy
:
    public option
{
    // Private data

        //- Name of velocity field; default = U
        word UName_;


public:

    //- Runtime type information
    TypeName("buoyancyEnergy");


    // Constructors

        //- Construct from explicit source name and mesh
        buoyancyEnergy
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        buoyancyEnergy(const buoyancyEnergy&) = delete;


    // Member Functions

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames) override;


        // Evaluate

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            ) override;


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict) override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const buoyancyEnergy&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
