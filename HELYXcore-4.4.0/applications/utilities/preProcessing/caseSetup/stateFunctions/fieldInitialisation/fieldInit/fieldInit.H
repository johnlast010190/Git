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
    (c) 2017-2024 Engys Ltd.

Class
    Foam::fieldInit

Description

SourceFiles
    fieldInit.C
    newPhaseChangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef fieldInit_H
#define fieldInit_H

#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "fvMesh/fvMesh.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class fieldInit Declaration
\*---------------------------------------------------------------------------*/

class fieldInit
{
protected:

    // Protected data

        //const reference to mesh
        const fvMesh& mesh_;

        const fvSolutionRegistry& localDb_;

        //field name
        const word fieldName_;

        //reference to field definition dictionary
        const dictionary fieldDict_;

        //init indicator
        bool initialized_;


public:

    //- Runtime type information
    TypeName("fieldInit");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fieldInit,
            initMethod,
            (
                const fvMesh& mesh,
                const fvSolutionRegistry& localDb,
                const dictionary& fieldDict,
                const word& fN
            ),
            (mesh, localDb, fieldDict, fN)
        );


    // Selectors

        //- Return a reference to the selected fieldInit model
        static autoPtr<fieldInit> New
        (
            const fvMesh& mesh,
            const fvSolutionRegistry& localDb,
            const dictionary& fieldDict,
            const word& fN
        );


    // Constructors

        //- Construct from components
        fieldInit
        (
            const fvMesh& mesh,
            const fvSolutionRegistry& localDb,
            const dictionary& fieldDict,
            const word& fN
        );

        //- Disallow default bitwise copy construction
        fieldInit(const fieldInit&) = delete;


    //- Destructor
    virtual ~fieldInit()
    {}


    // Member Functions

        // Access

            //- Return const access to mesh
            inline const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return const reference to the local database
            inline const objectRegistry& localDb() const
            {
                return localDb_.registry();
            }

            //- Return const reference to the local database
            inline const fvSolutionRegistry& solReg() const
            {
                return localDb_;
            }

            //- Return const access to field name
            inline word name() const
            {
                return fieldName_;
            }

            //- Return const access to field dictionary
            inline const dictionary& fieldDict() const
            {
                return fieldDict_;
            }

            //- Const access to modification switch
            inline bool initialised() const
            {
                return initialized_;
            }

            //- Non-const access to modification switch
            inline bool& initialised()
            {
                return initialized_;
            }

            //- Const access to initialisation information
            const dictionary& initDict() const;

        //- Print initialisation message
        void initMsg() const;

        //- Correct the fieldInit model
        virtual void correct() = 0;

        //- Init coupled boundaries to zero gradient
        virtual void initCoupledBoundaries() const;

        //- Init mixed and calculated type BCs from internal field
        virtual void initTurbBoundaries() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fieldInit&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
