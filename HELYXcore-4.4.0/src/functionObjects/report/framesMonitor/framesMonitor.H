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
    (c) 2023-2024 Engys Ltd.

Class
    Foam::functionObjects::framesMonitor

Group
    grpReportFunctionObjects

Description
    Prints state of the reference frame.

Usage
    \table
        Property     | Description             | Required    | Default value
        referenceFrame| name of the frame      | yes         | none
    \endtable

    Example of the boundary condition specification:
    \verbatim
    framesMonitor
    {
        type            framesMonitor;
        libs            ("libreport.so");
        writeControl    writeTime;
        referenceFrame rotorFrame;
    }
    \endverbatim

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    framesMonitor.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_framesMonitor_H
#define functionObjects_framesMonitor_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "db/Time/Time.H"
#include "cfdTools/general/include/fvCFD.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "containers/Lists/PackedList/PackedList.H"
#include "materialModels/materialTables/materialTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class framesMonitor Declaration
\*---------------------------------------------------------------------------*/

class framesMonitor
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

   // Protected data

        //- Reference frame name
        const word frameName_;

        //- Reference frame pointer
        const coordinateFrame& coorFrame_;


public:

    //- Runtime type information
    TypeName("framesMonitor");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        framesMonitor
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        framesMonitor(const framesMonitor&) = delete;


    //- Destructor
    virtual ~framesMonitor(){};


    // Member Functions

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Read the parameters
        virtual bool read(const dictionary&);

        //- Write the object into the file
        virtual void writeFileHeader(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const framesMonitor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
