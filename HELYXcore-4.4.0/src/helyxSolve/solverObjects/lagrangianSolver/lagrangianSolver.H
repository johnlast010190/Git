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
    (c) 2025 Engys Ltd.

Class
    Foam::fv::lagrangianSolver

Group
    grpFvOptionsSources

Description
    Thermal and species source terms from combustion model.

SourceFiles
    lagrangianSolver.C

\*---------------------------------------------------------------------------*/

#ifndef lagrangianSolver_H
#define lagrangianSolver_H

#include "solverObjects/solverObject/solverObject.H"
#include "rhoThermo/rhoThermo.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "clouds/derived/basicReactingCloud/basicReactingCloud.H"
#include "SLGThermo/SLGThermo.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "clouds/derived/basicReactingMultiphaseCloud/basicReactingMultiphaseCloud.H"
#include "clouds/derived/basicThermoCloud/basicThermoCloud.H"
#include "coalCloud/coalCloud.H"
#include "clouds/derived/basicSprayCloud/basicSprayCloud.H"
#include "clouds/derived/basicKinematicCollidingCloud/basicKinematicCollidingCloud.H"
#include "clouds/derived/basicKinematicCloud/basicKinematicCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class lagrangianSolver Declaration
\*---------------------------------------------------------------------------*/

class lagrangianSolver
:
    public solverObject
{
    // Private Member Functions

        //- Pickup or construct slgThermo
        const SLGThermo& slgThermo();


protected:

    // Protected Member Data

        //- Lagrangian parcels
        autoPtr<basicReactingCloud> parcels_;

        //- Lagrangian multiphase parcels
        autoPtr<basicReactingMultiphaseCloud> multiphseParcels_;

        //- Lagrangian coal cloud
        autoPtr<basicThermoCloud> basicThermoCloud_;

        //- Lagrangian multiphase parcels
        autoPtr<coalCloud> coalCloud_;

        //- Lagrangian basic spray cloud
        autoPtr<basicSprayCloud> basicSprayCloud_;

        //- The kinematic colliding cloud
        autoPtr<basicKinematicCollidingCloud> basicKinematicCollidingCloud_;

        //- The kinematic cloud
        autoPtr<basicKinematicCloud> basicKinematicCloud_;

        //- Pointer to solution control
        solutionControl* solnControlPtr_;

        //- Pointer to the fluid thermo
        rhoThermo *thermoPtr_;

        //- Pointer to gravitational acceleration
        autoPtr<uniformDimensionedVectorField> g_;

        //- Name of the lagrangian dict to pick up
        word cloudDictName_;


public:

    //- Runtime type information
    TypeName("lagrangianSolver");


    // Constructors

        //- Construct from components
        lagrangianSolver
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construct
        lagrangianSolver(lagrangianSolver&) = delete;


    //- Destructor
    virtual ~lagrangianSolver()
    {}

    //- Lagrangian cloud types
    enum lagrangianTypes
    {
        BasicReactingCloud,
        BasicReactingMultiphaseCloud,
        BasicThermoCloud,
        CoalCloud,
        BasicSprayCloud,
        BasicKinematicCollidingCloud,
        BasicKinematicCloud
    };

    //- Names of lagrangian types
    static const NamedEnum<lagrangianTypes, 7> lagrangianTypesNames_;

    //- Lagrangian type
    lagrangianTypes lagrangianType_;


    // Solver Functions

        virtual bool initialise() override;

        virtual void read(const dictionary& dict) override;

        virtual void getSolveGraph
        (
            wordList& solveNames,
            HashTable<wordList>& requiredDependencies,
            HashTable<wordList>& optionalDependencies,
            HashTable<wordList>& correctorMembers
        ) override;


        void beginIteration
        (
            const label corrector,
            const word& correctorName,
            const bool finalIter
        ) override;

        //- Add contribution to compressible equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldI
        ) override;

        //- Add contribution to equation
        void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldI
        ) override;

        //- Add contribution to compressible equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldI
        ) override;

        virtual void correct(const word&, const word&) override;


    // Source term functions

        virtual void getSourceGraph
        (
            wordList& fields,
            HashTable<wordList>& sourceDependencies
        ) override;

        scalar getMaxTimeStep() override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(lagrangianSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
