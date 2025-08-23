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
    (c) 2017 OpenFOAM Foundation
    (c) 2020-2024 Engys Ltd.

Class
    Foam::fv::radiationSolver

Description
    Calculates and applies the radiation source to the energy equation.

SourceFiles
    radiationSolver.C

\*---------------------------------------------------------------------------*/

#ifndef radiationSolver_H
#define radiationSolver_H

#include "solverObjects/solverObject/solverObject.H"
#include "solverObjects/solverOption/SolverOption.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "radiationModels/radiationModel/radiationModel.H"
#include "cfdTools/general/fvOptions/fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class transportModel;

namespace fv
{

/*---------------------------------------------------------------------------*\
               Class radiationSolver Declaration
\*---------------------------------------------------------------------------*/

class radiationSolver
:
    public solverObject
{
    // Private data

        //- The radiation model pointer
        Foam::radiationModel* radiationPtr_;

        //- Thermo model
        basicThermo* thermoPtr_;

        //- Transport model
        const transportModel* transportPtr_;

        //- Temperature field
        const volScalarField* TPtr_;

        //- For incompressible case
        dimensionedScalar rhoCpRef_;

        //- Constant part of the fvDOM qr field, saved at computation time
        //  Only the boundary is used, but need a GeometricField to get mapping
        //  calls
        autoPtr<volScalarField> qrConst_;


protected:

    void updateQrConst();


public:

    //- Runtime type information
    TypeName("radiationSolver");


    // Constructors

        //- Construct from components
        radiationSolver
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- Disallow copy construction
        radiationSolver(radiationSolver&) = delete;


    // Member Functions

        virtual bool initialise() override;

        virtual void getSolveGraph
        (
            wordList& solveNames,
            HashTable<wordList>& requiredDependencies,
            HashTable<wordList>& optionalDependencies,
            HashTable<wordList>& correctorMembers
        ) override;

        virtual void correct(const word&, const word&) override;

        virtual void getSourceGraph
        (
            wordList& fields,
            HashTable<wordList>& sourceDependencies
        ) override;

        //- Add contribution to temperature equation
        virtual void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldI
        ) override;

        //- Add contribution to energy equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldI
        ) override;

        virtual void getBoundarySourceGraph
        (
            HashTable<labelList>& fieldPatchIDs,
            HashTable<wordList>& boundarySourceDependencies
        ) override;


        //- Add boundary source and derivative
        virtual void addBoundarySource
        (
            const word& fieldName,
            const label patchID,
            const scalarField& pf,
            scalarField& f,
            scalarField& df
        ) override;


        //- Access to radiation model
        Foam::radiationModel& radiationModel()
        {
            if (!radiationPtr_)
            {
                FatalErrorInFunction
                    << "No radiation model set."
                    << exit(FatalError);
            }
            return *radiationPtr_;
        }


    // Member operators

        //- Disallow default bitwise assignment
        void operator=(radiationSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
