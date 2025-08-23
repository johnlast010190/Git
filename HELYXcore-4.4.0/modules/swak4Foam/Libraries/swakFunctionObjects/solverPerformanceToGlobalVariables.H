/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) ICE Stroemungsfoschungs GmbH
    (c) 2024 Engys Ltd.

Class
    solverPerformanceToGlobalVariables

Description
    Function object that gets the solverPerformances of specified
    solvers and stores them in global variables

SourceFiles
    solverPerformanceToGlobalVariables.C

Contributors/Copyright:
    2011, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2013 Bruno Santos <wyldckat@gmail.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#ifndef solverPerformanceToGlobalVariables_H
#define solverPerformanceToGlobalVariables_H

#include "fields/Fields/primitiveFieldsFwd.H"
#include "meshes/primitiveShapes/point/pointFieldFwd.H"

#include "include/swak.H"

namespace Foam
{
// Forward declaration of classes
class dictionary;
class objectRegistry;
class polyTopoChangeMap;
class polyMeshMap;
class polyMesh;

class solverPerformanceToGlobalVariables
{
    // Private Data

        const objectRegistry& obr_;

        //- Names of the Fields for which the values are added
        wordList fieldNames_;

        //- The global namespace the variables should go to
        word toGlobalNamespace_;

        //- Whether or not to reset the variable
        bool noReset_;


protected:

    // Protected Member Functions

        void addFieldToData(const word &name);

        void setValue(const word &name,scalar value);

        void setValue(const word &name,const scalarField &value);

        void executeAndWriteToGlobal();


public:

    //- Runtime type information
    TypeName("solverPerformanceToGlobalVariables");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        solverPerformanceToGlobalVariables
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );

        //- Disallow default bitwise copy construction
        solverPerformanceToGlobalVariables
        (
            const solverPerformanceToGlobalVariables&
        ) = delete;


    //- Destructor
    virtual ~solverPerformanceToGlobalVariables();


    // Member Functions

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Read the forces data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();
        virtual void setPostProcess(const bool&){};

        //- Execute at the final time-loop, currently does nothing
        virtual void end();


        //- Write the forces
#ifdef FOAM_IOFILTER_WRITE_NEEDS_BOOL
        virtual bool write();
#else
        virtual void write();
#endif


        //- Update for changes of mesh
        virtual void topoChange(const polyTopoChangeMap&)
        {}

        //- Update for changes of mesh
        virtual void mapMesh(const polyMeshMap&)
        {}

        //- Update for changes of mesh
#ifdef FOAM_MOVEPOINTS_GETS_POLYMESH
        virtual void movePoints(const polyMesh&)
#else
        virtual void movePoints(const pointField&)
#endif
        {}

        //- Clear all data
        void clearData();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solverPerformanceToGlobalVariables&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
