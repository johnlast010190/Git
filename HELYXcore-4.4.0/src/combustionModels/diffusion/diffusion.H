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
    (c) 2012-2022 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::combustionModels::diffusion

Group
    grpCombustionModels

Description
    Simple diffusion-based combustion model based on the principle mixed is
    burnt. Additional parameter C is used to distribute the heat release rate
    in time.

SourceFiles
    diffusion.C

\*---------------------------------------------------------------------------*/

#ifndef diffusion_H
#define diffusion_H

#include "singleStepCombustion/singleStepCombustion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                         Class diffusion Declaration
\*---------------------------------------------------------------------------*/

class diffusion
:
    public singleStepCombustion
{
    // Private data

        //- Model constant
        scalar C_;

        //- Name of oxidant - default is "O2"
        word oxidantName_;


public:

    //- Runtime type information
    TypeName("diffusion");


    // Constructors

        //- Construct from components
        diffusion
        (
            const word& modelType,
            const fluidMulticomponentThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );

        //- Disallow default bitwise copy construction
        diffusion(const diffusion&) = delete;


    //- Destructor
    virtual ~diffusion();


    // Member Functions

        //- Correct combustion rate
        virtual void correct();

        //- Update properties
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const diffusion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
