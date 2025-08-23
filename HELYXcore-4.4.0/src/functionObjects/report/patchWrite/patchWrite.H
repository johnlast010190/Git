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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2010-2024 Engys Ltd.

Class
    Foam::functionObjects::patchWrite

Group
    grpReportFunctionObjects

Description
    Write field value on selcted patches

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    patchWrite.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_patchWrite_H
#define functionObjects_patchWrite_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "fields/volFields/volFields.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class patchWrite Declaration
\*---------------------------------------------------------------------------*/

class patchWrite
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- input
        wordList writePatchNames_;
        labelList writePatchIDs_;
        wordList fieldNames_;

        Switch nearCellValue_;


public:

    //- Runtime type information
    TypeName("patchWrite");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        patchWrite
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        patchWrite(const patchWrite&) = delete;


    //- Destructor
    virtual ~patchWrite();


    // Member Functions

        //- Update solver
        void calculate();

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Read the parameters
        virtual bool read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const patchWrite&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
