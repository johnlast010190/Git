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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

Class
    Foam::combustionModels::eddyDissipationModel

Group
    grpCombustionModels

Description
    Simple infinitely fast chemistry combustion model based on the principle
    mixed is burnt. Additional parameter C is used to distribute the heat
    release rate.in time

SourceFiles
    eddyDissipationModel.C

\*---------------------------------------------------------------------------*/

#ifndef eddyDissipationModel_H
#define eddyDissipationModel_H

#include "singleStepCombustion/singleStepCombustion.H"
#include "LES/LESModel/LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class eddyDissipationModel Declaration
\*---------------------------------------------------------------------------*/

class eddyDissipationModel
:
    public singleStepCombustion
{
    // Private data

        //- Model constant
        scalar C_;

        scalar Cd_;

        scalar Cstiff_;

        //- Product volume - used by alphat BC

        volScalarField PV_;


    // Private Member Functions

        //- Return the reciprocal of the turbulent mixing time scale
        tmp<volScalarField> rtTurb() const;

        //- Return the reciprocal of the diffusion time scale
        tmp<volScalarField> rtDiff() const;

        //- Calculate the product volume
        void calcPV();


public:

    //- Runtime type information
    TypeName("eddyDissipation");


    // Constructors

        //- Construct from components
        eddyDissipationModel
        (
            const word& modelType,
            const fluidMulticomponentThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );

        //- Disallow default bitwise copy construction
        eddyDissipationModel(const eddyDissipationModel&) = delete;


    //- Destructor
    virtual ~eddyDissipationModel();


    // Member Functions

        //- Correct combustion rate
        virtual void correct();

        //- Update properties
        virtual bool read();


    // Member Operators

        //- No copy assignment
        void operator=(const eddyDissipationModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
