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
    (c) 2016 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::writeCellVolumes

Group
    grpFieldFunctionObjects

Description
    Writes the cell-volumes volScalarField.

Usage
    Example of function object specification:
    \verbatim
    writeCellVolumes
    {
        type        writeCellVolumes;
        libs        ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property  | Description                 | Required    | Default value
        type      | type name: writeCellVolumes | yes       |
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    writeCellVolumes.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writeCellVolumes_H
#define functionObjects_writeCellVolumes_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class writeCellVolumes Declaration
\*---------------------------------------------------------------------------*/

class writeCellVolumes
:
    public fvMeshFunctionObject
{
public:

    //- Runtime type information
    TypeName("writeCellVolumes");


    // Constructors

        //- Construct from Time and dictionary
        writeCellVolumes
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        writeCellVolumes(const writeCellVolumes&) = delete;


    //- Destructor
    virtual ~writeCellVolumes();


    // Member Functions

        //- Read the cell-volume data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Write the cell-volume fields
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const writeCellVolumes&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
