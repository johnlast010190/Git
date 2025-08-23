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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::processorField

Group
    grpFieldFunctionObjects

Description
    Writes a scalar field whose value is the local processor ID.  The output
    field name is 'processorID'.

Usage
    Example of function object specification:
    \verbatim
    processorField1
    {
        type        processorField;
        libs        ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | type name: processorField | yes       |
        log          | Log to standard output  | no          | yes
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    processorField.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_processorField_H
#define functionObjects_processorField_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class processorField Declaration
\*---------------------------------------------------------------------------*/

class processorField
:
    public fvMeshFunctionObject
{

        //- Result name
        word resultName_;

        //- Switch to send output to Info as well as to file
        Switch log_;


public:

    //- Runtime type information
    TypeName("processorField");


    // Constructors

        //- Construct from Time and dictionary
        processorField
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        processorField(const processorField&) = delete;


    //- Destructor
    virtual ~processorField();


    // Member Functions

        //- Read the input data
        virtual bool read(const dictionary&);

        //- Calculate the processorID field
        virtual bool execute();

        //- Write the processorID field
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const processorField&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
