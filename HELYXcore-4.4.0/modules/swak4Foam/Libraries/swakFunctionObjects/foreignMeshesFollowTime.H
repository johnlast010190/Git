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
    foreignMeshesFollowTime

Description
    Function object that makes a number of foreign meshes "follow" the
    current time

SourceFiles
    foreignMeshesFollowTime.C

Contributors/Copyright:
    2012-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2013 Bruno Santos <wyldckat@gmail.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef foreignMeshesFollowTime_H
#define foreignMeshesFollowTime_H

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

class foreignMeshesFollowTime
{
    // Private Data

        //- Names of the meshes that should be followed
        wordList followingMeshes_;

        //- Do the actual work
        void followCurrentTime();

        //- Object registry
        const objectRegistry &obr_;


public:

    //- Runtime type information
    TypeName("foreignMeshesFollowTime");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        foreignMeshesFollowTime
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );

        //- Disallow default bitwise copy construction
        foreignMeshesFollowTime(const foreignMeshesFollowTime&) = delete;


    //- Destructor
    virtual ~foreignMeshesFollowTime();


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
        void operator=(const foreignMeshesFollowTime&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
