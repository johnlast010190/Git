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
    Foam::fv::surfaceFilmSolver

Group
    grpFvOptionsSources

Description
    Thermal and species source terms from combustion model.

SourceFiles
    surfaceFilmSolver.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceFilmSolver_H
#define surfaceFilmSolver_H

#include "solverObjects/solverObject/solverObject.H"
#include "surfaceFilmModels/surfaceFilmModel/surfaceFilmModel.H"
#include "rhoThermo/rhoThermo.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class surfaceFilmSolver Declaration
\*---------------------------------------------------------------------------*/

class surfaceFilmSolver
:
    public solverObject
{
protected:


    // Protected Member Data

        //- Surface film model
        autoPtr<regionModels::surfaceFilmModel> surfaceFilm_;

        //- Pointer to solution control
        solutionControl* solnControlPtr_;

        //- Pointer to the fluid thermo
        rhoThermo *thermoPtr_;

        //- Pointer to gravitational acceleration
        autoPtr<uniformDimensionedVectorField> g_;


public:

    //- Runtime type information
    TypeName("surfaceFilmSolver");


    // Constructors

        //- Construct from components
        surfaceFilmSolver
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construct
        surfaceFilmSolver(surfaceFilmSolver&) = delete;


    //- Destructor
    virtual ~surfaceFilmSolver()
    {}


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

        virtual void correct(const word&, const word&) override;


    // Source term functions

        virtual void getSourceGraph
        (
            wordList& fields,
            HashTable<wordList>& sourceDependencies
        ) override;

        //- Add contribution to compressible equation
        void addSup(fvMatrix<scalar>& eqn, const label fieldI) override;

        //- Add contribution to compressible equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldI
        ) override;

        scalar getMaxTimeStep() override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(surfaceFilmSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
