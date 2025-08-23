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
    (c) 2019 Engys Ltd.

Class
    Foam::fv::SolverFunction

Description
    A function object template class to instantiate solverObjects as function
    objects.

SourceFiles
    SolverFunction.C

\*---------------------------------------------------------------------------*/

#ifndef SolverFunction_H
#define SolverFunction_H

#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatricesFwd.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "db/dictionary/dictionary.H"
#include "primitives/bools/Switch/Switch.H"
#include "solverObjects/solverObject/solverObject.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"
#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/IOstreams/Fstreams/OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class SolverFunction Declaration
\*---------------------------------------------------------------------------*/

template<class SolverObject>
class SolverFunction
:
    public fvMeshFunctionObject
{

protected:

    // Protected data

        SolverObject solverObject_;

        bool postProcess_;

public:

    //- Runtime type information
    TypeName("SolverFunction");


    // Constructor

        //- Construct from components
        SolverFunction
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool loadFromFiles=false
        );

    // Destructor

        virtual ~SolverFunction();

    // Member Functions

        //- Read source dictionary
        virtual bool read(const dictionary& dict);

        //- Execute
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Execute at the final time loop
        virtual bool end();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSolverFunction(SolverObject)                                       \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Foam::functionObjects::SolverFunction<Foam::fv::SolverObject>,                      \
        #SolverObject,                                                         \
        0                                                                      \
    );                                                                         \
                                                                               \
    Foam::functionObject::adddictionaryConstructorToTable                      \
        <Foam::functionObjects::SolverFunction<Foam::fv::SolverObject>>        \
        addSolverFunction##SolverObject##dictionary##ConstructorTofunctionObjectTable_


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "SolverFunction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
