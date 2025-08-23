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
    Foam::combustionModels::eddyDissipationModelBase

Description
    Standard Eddy Dissipation Model based on the assumption that the
    reaction rates are controlled by the turbulece time length scale.
    The chemical reaction rate is governed by the large-eddy mixing scale,
    (ie k/epsilon). Combustion takes place whenever the turbulence is present.


SourceFiles
    eddyDissipationModelBase.C

\*---------------------------------------------------------------------------*/

#ifndef eddyDissipationModelBase_H
#define eddyDissipationModelBase_H

#include "singleStepCombustion/singleStepCombustion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class eddyDissipationModelBase Declaration
\*---------------------------------------------------------------------------*/

class eddyDissipationModelBase
:
    public singleStepCombustion
{
    // Private data

        //- EDM constant
        scalar CEDC_;


public:

    // Constructors

        //- Construct from components
        eddyDissipationModelBase
        (
            const word& modelType,
            const fluidMulticomponentThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );

        //- Disallow default bitwise copy construction
        eddyDissipationModelBase(const eddyDissipationModelBase&) = delete;


    //- Destructor
    virtual ~eddyDissipationModelBase();


    // Member Functions

        //- Return the reciprocal of the turbulent mixing time scale
        tmp<volScalarField> rtTurb() const;

        //- Correct combustion rate
        void correct();

        //- Calculate time scale
        virtual tmp<volScalarField> timeScale() = 0;

        //- Update properties
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const eddyDissipationModelBase&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "eddyDissipationModelBase/eddyDissipationModelBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
