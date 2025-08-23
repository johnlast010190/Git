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
    (c) 2015 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::radiationModels::transmissivityModels::constantTransmissivity

Description
    Constant radiation transmissivity coefficient

SourceFiles
    constantTransmissivity.C

\*---------------------------------------------------------------------------*/

#ifndef radiationConstantTransmissivity_H
#define radiationConstantTransmissivity_H

#include "transmissivityModels/transmissivityModel/transmissivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace transmissivityModels
{

/*---------------------------------------------------------------------------*\
                   Class constantTransmissivity Declaration
\*---------------------------------------------------------------------------*/

class constantTransmissivity
:
    public transmissivityModel
{

    // Private data

        //- Coefficients dictionary
        dictionary coeffsDict_;

        //- Transmissivity coefficient / []
        scalar tau_;


public:

    //- Runtime type information
    TypeName("constantTransmissivity");


    // Constructors

        //- Construct from components
        constantTransmissivity(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~constantTransmissivity();


    // Member Functions

        //- Return scatter coefficient
        tmp<volScalarField> tauEff(const label bandI = 0) const;


        //- Is Grey
        inline bool isGrey() const
        {
            return true;
        }


        //- Number of bands
        inline label nBands() const
        {
            return 1;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace transmissivityModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
