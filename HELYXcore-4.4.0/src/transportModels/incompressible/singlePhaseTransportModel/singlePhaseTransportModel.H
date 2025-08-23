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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2010-2024 Engys Ltd.

Class
    Foam::singlePhaseTransportModel

Description
    A simple single-phase transport model based on viscosityModel.

    Used by the incompressible single-phase solvers like simpleFoam,
    turbFoam etc.

SourceFiles
    singlePhaseTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef singlePhaseTransportModel_H
#define singlePhaseTransportModel_H

#include "transportModel/transportModel.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "physicalProperties/physicalProperties.H"
#include "memory/autoPtr/autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class viscosityModel;

/*---------------------------------------------------------------------------*\
                Class singlePhaseTransportModel Declaration
\*---------------------------------------------------------------------------*/

class singlePhaseTransportModel
:
    public IOdictionary,
    public transportModel
{
    // Private Data

        autoPtr<viscosityModel> viscosityModelPtr_;
        autoPtr<physicalProperties> physicalPropertiesPtr_;
        mutable autoPtr<volScalarField> Cp_;


    // Private Member Functions

        //- return const access to physicalProperties
        const physicalProperties& properties() const;

        //- eturn const access to viscosityModel
        const viscosityModel& viscosity() const;


public:

    //- Runtime type information
    TypeName("singlePhaseTransportModel");


    // Constructors

        //- Construct from components
        singlePhaseTransportModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );

        //- Disallow default bitwise copy construction
        singlePhaseTransportModel(const singlePhaseTransportModel&) = delete;


    //- Destructor
    virtual ~singlePhaseTransportModel();


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const;

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const;

        //- Return constant density
        tmp<volScalarField> rho() const;

        //- Constant heat capacity
        const volScalarField& Cp() const;

        //- Constant thermal conductivity
        tmp<volScalarField> lambda() const;

        //- Constant turbulent Prandtl number
        tmp<volScalarField> Prt() const;

        //- Laminar Prandtl number
        tmp<volScalarField> Pr() const;

        //- Laminar thermal diffusivity
        tmp<volScalarField> alphaLam() const;

        //- Correct the laminar viscosity
        virtual void correct();

        //- Read transportProperties dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const singlePhaseTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
