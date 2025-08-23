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
    (c) 2016-2024 Engys Ltd.
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2017 OpenCFD Ltd

Class
    Foam::solutionControl

Description
    Base class for solution control classes

\*---------------------------------------------------------------------------*/

#ifndef solutionControl_H
#define solutionControl_H

#include "fvMesh/fvMesh.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solutionControl Declaration
\*---------------------------------------------------------------------------*/

class solutionControl
:
    public regIOobject
{
public:

    struct fieldData
    {
        wordRe name;
        scalar absTol;
        scalar relTol;
        scalar initialResidual;
    };


    // Static member functions

        template<class Type>
        static void maxTypeResidual
        (
            const objectRegistry& obr,
            const entry& solverPerfDictEntry,
            scalar& firstRes,
            scalar& lastRes
        );

        template<class Type>
        static void maxTypeBlockResidual
        (
            const objectRegistry& obr,
            const entry& solverPerfDictEntry,
            scalar& firstRes,
            scalar& lastRes
        );

        static scalar maxResidual
        (
            const objectRegistry& obr,
            const entry& solverPerfDictEntry,
            scalar& lastRes
        );

        //- Return index of field in residualControl_ if present
        virtual label applyToField
        (
            const word& fieldName,
            const bool useRegEx = true,
            const word& correctorName = word::null
        ) const;



protected:

    // Protected data

        //- Reference to the mesh
        fvMesh& mesh_;

        //- Reference to the mesh database
        const objectRegistry& obr_;

        //- List of residual data per field, for each corrector
        //  Backward-compatible data is keyed with word::null
        HashTable<List<fieldData>> residualControl_;

        //- Minimum iterations for each corrector
        HashTable<label> minIters_;

        //- Requested iterations for each corrector
        HashTable<label> nIters_;

        //- The dictionary name, e.g. SIMPLE, PIMPLE
        const word algorithmName_;


        // Solution controls

            //- Maximum number of non-orthogonal correctors
            label nNonOrthCorr_;

            //- Flag to indicate to solve for momentum
            bool momentumPredictor_;

            //- Flag to indicate to solve using transonic algorithm
            bool transonic_;

            //- Flag to indicate to relax pressure using the
            //  "consistent" approach of SIMPLEC
            bool consistent_;

            //- Flag for Choi's modified momentum interpolation method
            bool modifiedMomentumInterp_;

            //- Flag to indicate that the flow system of equations should not
            //  be evolved
            bool frozenFlow_;

            //- Flag to indicate that an adjoint optimization loop is active
            bool adjointOptimization_;

        // Evolution

            //- Current corrector loop index
            label corr_;

            //- Current non-orthogonal corrector loop index
            label corrNonOrtho_;


    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        virtual void readControls(const bool absTolOnly);

        //- Read residual controls for a certain corrector from fvSolution
        void readResidualControls
        (
            const word& correctorName,
            const dictionary& residualDict,
            const bool absTolOnly
        );

        //- Return true if equation variables should be stored
        virtual bool storeVars() const;

        //- Store previous iteration field for vol<Type>Fields
        template<class Type>
        void storePrevIter() const;

        //- Set the firstIteration flag on the mesh data dictionary
        virtual void setFirstIterFlag
        (
            const bool check = true,
            const bool force = false
        );


public:


    // Static Data Members

        //- Run-time type information
        TypeName("solutionControl");

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            solutionControl,
            dictionary,
            (fvMesh& mesh, const objectRegistry& obr),
            (mesh, obr)
        );


    // Constructors

        //- Construct from registry
        solutionControl(fvMesh& mesh, const word& algorithmName);

        //- Construct from registry
        solutionControl
        (
            fvMesh& mesh,
            const objectRegistry& obr,
            const word& algorithmName
        );

        //- Disallow default bitwise copy construction
        solutionControl(const solutionControl&) = delete;


    // Selectors

        //- Run-time selector
        static autoPtr<solutionControl> New
        (
            fvMesh& mesh, const objectRegistry& obr
        );

        static autoPtr<solutionControl> NewSIMPLE
        (
            fvMesh& mesh, const objectRegistry& obr
        );

        //- Look up or construct if not found
        static solutionControl& lookupOrCreate
        (
            fvMesh& mesh, const objectRegistry& obr
        );

        //- Look up or construct a specific derived type
        template<class Type, class ...Args>
        static Type& lookupOrCreateType
        (
            fvMesh& mesh, const objectRegistry& obr, Args&&... args
        );

    //- Destructor
    virtual ~solutionControl();


    // Member Functions

        // Access

            //- Return the solution dictionary
            inline const dictionary& dict() const;

            //- Current corrector loop index
            inline label corr() const;

            //- Set current corrector loop index
            inline void setCorr(const label corr);

            //- Current non-orthogonal corrector index
            inline label corrNonOrtho() const;


        // Solution control

            //- Maximum number of outer correctors per timestep
            virtual label nOuterCorr() const
            {
                return 1;
            }

            //- Maximum number of non-orthogonal correctors
            inline label nNonOrthCorr() const;

            //- Helper function to identify final non-orthogonal iteration
            inline bool finalNonOrthogonalIter() const;

            //- Flag to indicate to solve for momentum
            inline bool momentumPredictor() const;

            //- Flag to indicate to solve using transonic algorithm
            inline bool transonic() const;

            //- Flag to indicate to relax pressure using the
            //  "consistent" approach of SIMPLEC
            inline bool consistent() const;

            //- Flag to indicate that the flow system of equations should not
            //  be evolved
            inline bool frozenFlow() const;

            //- Flag for Choi's modified momentum interpolation method
            inline bool modifiedMomentumInterp() const;

            inline bool adjointOptimization() const;

            void setAdjointOptimization();

            const List<fieldData>& residualControl
            (
                const word& correctorName = word::null
            ) const
            {
                return residualControl_[correctorName];
            }

            List<fieldData>& residualControl
            (
                const word& correctorName = word::null
            )
            {
                return residualControl_[correctorName];
            }

            scalar maxResidual
            (
                const entry& solverPerfDictEntry,
                scalar& lastRes
            ) const;

            //- Minium number of iterations requested for each corrector
            //  (this number of correctors will always be performed)
            label minIter(const word& correctorName) const
            {
                return minIters_.lookup(correctorName, 0);
            }

            //- Number of iterations requested for each corrector
            //  (less will be performed if converged; more if requested by
            //  solver objs)
            label nIter(const word& correctorName) const
            {
                return nIters_.lookup(correctorName, 1);
            }


        // Evolution

            //- Read controls from fvSolution dictionary
            virtual void readControls();

            //- Store previous iteration fields
            virtual void storePrevIterFields() const;

            //- Main control loop
            virtual bool loop() = 0;

            //- Non-orthogonal corrector loop
            inline bool correctNonOrthogonal();

            //- Return true if all convergence checks are satisfied
            virtual bool criteriaSatisfied() = 0;

        // Writing

            //- regIOobject write data function
            virtual bool writeData(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solutionControl&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cfdTools/general/solutionControl/solutionControl/solutionControlTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cfdTools/general/solutionControl/solutionControl/solutionControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
