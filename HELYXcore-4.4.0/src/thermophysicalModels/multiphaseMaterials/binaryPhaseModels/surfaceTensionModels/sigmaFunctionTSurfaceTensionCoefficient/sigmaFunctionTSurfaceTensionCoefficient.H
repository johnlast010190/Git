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
    (c) 2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd

Class
    Foam::surfaceTensionModels::sigmaFunctionTSurfaceTensionCoefficient

Description
    Surface tension model with sigma as a Function1 of Temperature.

SourceFiles
    sigmaFunctionTSurfaceTensionCoefficient.C

\*---------------------------------------------------------------------------*/

#ifndef sigmaFunctionTSurfaceTensionCoefficient_H
#define sigmaFunctionTSurfaceTensionCoefficient_H

#include "../surfaceTensionModel/surfaceTensionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class basicThermo;

namespace surfaceTensionModels
{

/*---------------------------------------------------------------------------*\
              Class sigmaFunctionTSurfaceTensionCoefficient Declaration
\*---------------------------------------------------------------------------*/

class sigmaFunctionTSurfaceTensionCoefficient
:
    public surfaceTensionModel
{
    // Private data

        //- Surface tension value
        autoPtr<Function1<scalar>> sigma_;

        const basicThermo& thermo_;

public:

    //- Runtime type information
    TypeName("sigmaFunctionT");


    // Constructors

        //- Construct from a dictionary and a phase pair
        sigmaFunctionTSurfaceTensionCoefficient
        (
            const dictionary& dict,
            const phasePairKey& pair,
            const fvMesh& mesh,
            const bool registerObject
        );


    //- Destructor
    virtual ~sigmaFunctionTSurfaceTensionCoefficient();


    // Member Functions

        //- Return surface tension field
        virtual tmp<volScalarField> sigma() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceTensionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
