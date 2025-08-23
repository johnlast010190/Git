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
    manipulateField

Description
    Function object that manipulates a field

SourceFiles
    manipulateField.C

Contributors/Copyright:
    2010-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2013 Bruno Santos <wyldckat@gmail.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef manipulateField_H
#define manipulateField_H

#include "fields/Fields/primitiveFieldsFwd.H"
#include "fields/volFields/volFieldsFwd.H"
#include "meshes/primitiveShapes/point/pointFieldFwd.H"
#include "primitives/bools/Switch/Switch.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/dictionary/dictionary.H"

#include "include/swak.H"

#include "helpers/exprString.H"

namespace Foam
{
// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyTopoChangeMap;
class polyMeshMap;
class polyMesh;

class FieldValueExpressionDriver;

class manipulateField
{
protected:

    // Protected data

        //- Are we working?
        bool active_;

        //- Re-write manipulated field
        bool writeManipulated_;

        //- Name of the field
        word name_;

        //- Expression
        exprString expression_;

        //- Expression that says which cells are to be changed
        exprString maskExpression_;

        const objectRegistry& obr_;

        dictionary dict_;

        autoPtr<FieldValueExpressionDriver> driver_;


    // Protected Member Functions

        template<class TData, class TMask>
        void manipulate
        (
            const TData& data,
            const TMask& mask,
            const word entity = "cells"
        );

        template<class TData, class TMask>
        void manipulateSurface(const TData& data, const TMask& mask);


public:

    //- Runtime type information
    TypeName("manipulateField");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        manipulateField
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );

        //- Disallow default bitwise copy construction
        manipulateField(const manipulateField&) = delete;


    //- Destructor
    virtual ~manipulateField();


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
        void operator=(const manipulateField&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
