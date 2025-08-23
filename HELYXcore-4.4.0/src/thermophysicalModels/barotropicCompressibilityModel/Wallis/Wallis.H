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
    (c) 2011-2019 OpenFOAM Foundation

Class
    Foam::compressibilityModels::Wallis

Description
    Wallis compressibility model.

SourceFiles
    Wallis.C

\*---------------------------------------------------------------------------*/

#ifndef Wallis_H
#define Wallis_H

#include "barotropicCompressibilityModel/barotropicCompressibilityModel.H"
#include "dimensionedTypes/dimensionedScalar/dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressibilityModels
{

/*---------------------------------------------------------------------------*\
                           Class Wallis Declaration
\*---------------------------------------------------------------------------*/

class Wallis
:
    public barotropicCompressibilityModel
{
    // Private data

        dimensionedScalar pSat_;
        dimensionedScalar psiv_;
        dimensionedScalar psil_;

        dimensionedScalar rhovSat_;
        dimensionedScalar rholSat_;


public:

    //- Runtime type information
    TypeName("Wallis");


    // Constructors

        //- Construct from components
        Wallis
        (
            const dictionary& compressibilityProperties,
            const volScalarField& gamma,
            const word& psiName = "psi"
        );


    //- Destructor
    ~Wallis()
    {}


    // Member Functions

        //- Correct the Wallis compressibility
        void correct();

        //- Read transportProperties dictionary
        bool read(const dictionary& compressibilityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressibilityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
