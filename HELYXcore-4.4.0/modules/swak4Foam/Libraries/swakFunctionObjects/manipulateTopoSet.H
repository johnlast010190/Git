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
    manipulateTopoSet

Description
    Function object that manipulates a topoSet

SourceFiles
    manipulateTopoSet.C

Contributors/Copyright:
    2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef manipulateTopoSet_H
#define manipulateTopoSet_H

#include "primitives/bools/Switch/Switch.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/dictionary/dictionary.H"

#include "helpers/exprString.H"
#include "fields/Fields/scalarField/scalarField.H"
#include "meshes/primitiveShapes/point/pointField.H"

#include "include/swak.H"

namespace Foam
{
// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyTopoChangeMap;
class polyMeshMap;
class polyMesh;

class FieldValueExpressionDriver;

template<class TopoSetType>
class manipulateTopoSet
{
    // Private Member Functions

        void toMask(const scalarField &);

protected:

    // Protected data

        //- Are we working?
        bool active_;

        //- Re-write manipulated topoSet
        bool writeManipulated_;

        //- Should a non-existing cellSet be created
        bool createMissing_;

        //- Name of the topoSet
        word name_;

        //- Expression that says which are part of the topoSet
        exprString maskExpression_;

        const objectRegistry& obr_;

        dictionary dict_;

        //- Mask with the entities to add
        List<bool> mask_;

        autoPtr<FieldValueExpressionDriver> driver_;


    // Protected Member Functions

        //- This part is different for every topoSet
        void calcMask();


public:

    //- Runtime type information
    TypeName("manipulateTopoSet");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load topoSets from files
        manipulateTopoSet
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );

        //- Disallow default bitwise copy construction
        manipulateTopoSet(const manipulateTopoSet&) = delete;


    //- Destructor
    virtual ~manipulateTopoSet();


    // Member Functions

        //- Return name of the set of forces
        virtual const word& name() const
        {
            return name_;
        }

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


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const manipulateTopoSet&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
