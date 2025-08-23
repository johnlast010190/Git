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
    (c) 2016 OpenCFD Ltd
    (c) 2024-2025 Engys Ltd.

Class
    Foam::eddyDissipationDiffusionModel

Group
    grpCombustionModels

Description
    Eddy dissipation model based on the principle of mixed is burnt.

    The time scale is determined by the minimum of the turbulence
    and diffusion time scales.


SourceFiles
    eddyDissipationDiffusionModel.C

\*---------------------------------------------------------------------------*/

#ifndef eddyDissipationDiffusionModel_H
#define eddyDissipationDiffusionModel_H

#include "eddyDissipationModelBase/eddyDissipationModelBase.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                Class eddyDissipationDiffusionModel Declaration
\*---------------------------------------------------------------------------*/

class eddyDissipationDiffusionModel
:
    public eddyDissipationModelBase
{
    // Private data

        //- Diffussivity constant
        scalar Cd_;


public:

    //- Runtime type information
    TypeName("EDM");


    // Constructors

        //- Construct from components
        eddyDissipationDiffusionModel
        (
            const word& modelType,
            const fluidMulticomponentThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );

        //- Disallow default bitwise copy construction
        eddyDissipationDiffusionModel
        (
            const eddyDissipationDiffusionModel&
        ) = delete;


    //- Destructor
    virtual ~eddyDissipationDiffusionModel();


    // Member Functions

        //- Calculate time scale
        virtual tmp<volScalarField> timeScale();

         //- Return the reciprocal of the diffusion time scale
        tmp<volScalarField> rtDiff() const;

        //- Update properties
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const eddyDissipationDiffusionModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
