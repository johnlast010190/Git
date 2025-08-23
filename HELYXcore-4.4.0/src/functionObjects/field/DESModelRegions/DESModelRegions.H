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
    (c) 2016 OpenCFD Ltd.
    (c) 2013-2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::DESModelRegions

Group
    grpFieldFunctionObjects

Description
    This function object writes out an indicator field for DES turbulence
    calculations, that is:
    - 0 for RAS regions
    - 1 for LES regions

    The field is stored on the mesh database so that it can be retrieved and
    used for other applications.

Usage
    Example of function object specification to generate DES indicator field:
    \verbatim
    DESModelRegions1
    {
        type        DESModelRegions;
        libs        ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | type name: DESModelRegions| yes       |
        resultName   | Name of DES indicator field | no      | <function name>
        log          | log to standard output  | no          | yes
    \endtable

SourceFiles
    DESModelRegions.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_DESModelRegions_H
#define functionObjects_DESModelRegions_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "fields/volFields/volFieldsFwd.H"
#include "db/IOstreams/Fstreams/OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class DESModelRegions Declaration
\*---------------------------------------------------------------------------*/

class DESModelRegions
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        //- Result name
        word resultName_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("DESModelRegions");


    // Constructors

        //- Construct from Time and dictionary
        DESModelRegions
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        DESModelRegions(const DESModelRegions&) = delete;


    //- Destructor
    virtual ~DESModelRegions();


    // Member Functions

        //- Read the DESModelRegions data
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Calculate the DESModelRegions and write
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const DESModelRegions&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
