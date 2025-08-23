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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2024-2025 Engys Ltd.

Class
    Foam::combustionModel

Group
    grpcombustionModels

Description
    Base class for combustion models

SourceFiles
    combustionModel.C

\*---------------------------------------------------------------------------*/

#ifndef combustionModel_H
#define combustionModel_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "fluidMulticomponentThermo/fluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class combustionModel Declaration
\*---------------------------------------------------------------------------*/

class combustionModel
:
    public IOdictionary
{
    // Private Member Functions

        //- Construct the base IO object
        IOobject createIOobject
        (
            const fluidMulticomponentThermo& thermo,
            const word& combustionProperties
        ) const;


protected:

    // Protected data

        //- Reference to solution database
        const objectRegistry& obr_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the thermo
        const fluidMulticomponentThermo& thermo_;

        //- Reference to the turbulence model
        const compressibleTurbulenceModel& turb_;

        //- Dictionary of the model
        dictionary coeffs_;

        //- Model type
        const word modelType_;


public:

    //- Runtime type information
    TypeName("combustionModel");

    //- Default combustionProperties dictionary name
    static const word combustionPropertiesName;


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        combustionModel,
        dictionary,
        (
            const word& modelType,
            const fluidMulticomponentThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        ),
        (modelType, thermo, turb, combustionProperties)
    );


    // Constructors

        //- Construct from components
        combustionModel
        (
            const word& modelType,
            const fluidMulticomponentThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties = combustionPropertiesName
        );

        //- Disallow default bitwise copy construction
        combustionModel(const combustionModel&) = delete;


    // Selectors

        //- Generic New for each of the related chemistry model
        static autoPtr<combustionModel> New
        (
            const fluidMulticomponentThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties = combustionPropertiesName
        );


    //- Destructor
    virtual ~combustionModel();


    // Member Functions

        //- Return const access to the mesh database
        inline const fvMesh& mesh() const;

        //- Return const access to the thermo
        inline const fluidMulticomponentThermo& thermo() const;

        //- Return access to turbulence
        inline const compressibleTurbulenceModel& turbulence() const;

        //- Return const access to rho
        inline const volScalarField& rho() const;

        //- Return const access to phi
        inline tmp<surfaceScalarField> phi() const;

        //- Return const dictionary of the model
        inline const dictionary& coeffs() const;

        //- Correct combustion rate
        virtual void correct() = 0;

        //- Fuel consumption rate matrix, i.e. source term for fuel equation
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const = 0;

        //- Heat release rate [kg/m/s3]
        virtual tmp<volScalarField> Qdot() const = 0;

        //- Update properties from given dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const combustionModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "combustionModel/combustionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
