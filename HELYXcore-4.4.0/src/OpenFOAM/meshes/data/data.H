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
    Foam::data

Description
    Database for solution data, solver performance and other reduced data.

    fvMesh is derived from data so that all fields have access to the data from
    the mesh reference they hold.

SourceFiles
    data.C

\*---------------------------------------------------------------------------*/

#ifndef data_H
#define data_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "matrices/LduMatrix/LduMatrix/solverPerformance.H"
#include "matrices/blockLduMatrix/BlockLduSolvers/BlockLduSolver/BlockSolverPerformance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class data Declaration
\*---------------------------------------------------------------------------*/

class data
:
    public IOdictionary
{
    // Private data

        //- Previously used time-index, used for reset between iterations
        mutable DynamicList<Tuple2<word, label>> prevTimeIndices_;

        //- Previously used time-index, used for reset between iterations
        //  for the block residual dictionary
        mutable DynamicList<Tuple2<word, label>> prevTimeBlockIndices_;


public:

    //- Debug switch
    static int debug;


    // Constructors

        //- Construct for objectRegistry
        data(const objectRegistry& obr);

        //- Move constructor
        data(data&&) = default;

        //- Disallow default bitwise copy construction
        data(const data&) = delete;


    // Member Functions

        // Access

            //- Return the dictionary of solver performance data
            //  which includes initial and final residuals for convergence
            //  checking
            const dictionary& solverPerformanceDict() const;

            //- Add/set the solverPerformance entry for the named field
            template<class Type>
            void setSolverPerformance
            (
                const word& name,
                const SolverPerformance<Type>&
            ) const;

            //- Add/set the solverPerformance entry, using its fieldName
            template<class Type>
            void setSolverPerformance
            (
                const SolverPerformance<Type>&
            ) const;

            //- Return the dictionary of solver performance data
            //  which includes initial and final residuals for convergence
            //  checking
            const dictionary& blockSolverPerformanceDict() const;

            //- Add/set the solverPerformance entry for the named field
            template<class Type>
            void setBlockSolverPerformance
            (
                const word& name,
                const BlockSolverPerformance<Type>&
            ) const;

            //- Add/set the solverPerformance entry, using its fieldName
            template<class Type>
            void setBlockSolverPerformance
            (
                const BlockSolverPerformance<Type>&
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const data&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "meshes/data/dataTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
