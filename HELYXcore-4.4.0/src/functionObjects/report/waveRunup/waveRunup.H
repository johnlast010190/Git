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
    (c) 2016-2024 Engys Ltd.

Class
    Foam::functionObjects::waveRunup

Group
    grpReportFunctionObjects

Description
    Calculate height of wave against a patch. Actually works off the adjacent
    layer of cells in order to treat parallel boundaries correctly.

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    waveRunup.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_waveRunup_H
#define functionObjects_waveRunup_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "db/Time/Time.H"
#include "cfdTools/general/include/fvCFD.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "containers/Lists/PackedList/PackedList.H"
#include "convergenceTermination/convergenceTermination.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class waveRunup Declaration
\*---------------------------------------------------------------------------*/

class waveRunup
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private data

        //- Read from dictonary
        word alphaName_;
        wordList patchNames_;
        scalar contourLevel_;

        //- Result values
        scalar volHeight_;
        scalarList patchHeights_;


public:

    //- Runtime type information
    TypeName("waveRunup");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        waveRunup
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        waveRunup(const waveRunup&) = delete;


    //- Destructor
    virtual ~waveRunup();


    // Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os) const;

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
        void operator=(const waveRunup&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
