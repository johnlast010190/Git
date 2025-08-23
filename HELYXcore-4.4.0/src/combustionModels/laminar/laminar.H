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
    (c) 2013-2022 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::combustionModels::laminar

Group
    grpCombustionModels

Description
    Laminar combustion model.

SourceFiles
    laminar.C

\*---------------------------------------------------------------------------*/

#ifndef combustionModels_laminar_H
#define combustionModels_laminar_H

#include "memory/autoPtr/autoPtr.H"
#include "combustionModel/combustionModel.H"
#include "basicChemistryModel/basicChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                            Class laminar Declaration
\*---------------------------------------------------------------------------*/

class laminar
:
    public combustionModel
{
    // Private data

        //- Integrate reaction rate over the time-step
        //  using the selected ODE solver. Default true.
        bool integrateReactionRate_;

        //- Maximum integration time permitted. Default vGreat.
        scalar maxIntegrationTime_;

        //- Run chemistry rate calculation on every outer iteration. Default
        //  false. Has no effect if integrating.
        bool outerCorrect_;

        //- The time index of the last correction
        label timeIndex_;


protected:

    // Protected Data

        //- Pointer to chemistry model
        autoPtr<basicChemistryModel> chemistryPtr_;


public:

    //- Runtime type information
    TypeName("laminar");


    // Constructors

        //- Construct from components
        laminar
        (
            const word& modelType,
            const fluidMulticomponentThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );

        //- Disallow default bitwise copy construction
        laminar(const laminar&) = delete;


    //- Destructor
    virtual ~laminar();


    // Member Functions

        //- Correct combustion rate
        virtual void correct();

        //- Fuel consumption rate matrix.
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Heat release rate [kg/m/s3]
        virtual tmp<volScalarField> Qdot() const;

        //- Update properties from given dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const laminar&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
