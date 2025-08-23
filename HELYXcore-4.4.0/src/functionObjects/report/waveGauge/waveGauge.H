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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2016-2024 Engys Ltd.

Class
    Foam::waveGauge

Description
    Set of sets to sample.
    Call waveGauge.write() to sample&write files.

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    waveGauge.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_waveGauge_H
#define functionObjects_waveGauge_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "sampledSet/sampledSets/sampledSets.H"
#include "fields/volFields/volFieldsFwd.H"
#include "meshSearch/meshSearch.H"
#include "interpolation/interpolation/interpolation/interpolation.H"
#include "coordSet/coordSet.H"
#include "sampledSetWriters/writer.H"
#include "primitives/strings/lists/wordReList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class waveGauge Declaration
\*---------------------------------------------------------------------------*/

class waveGauge
:
    public sampledSets,
    public writeFile
{
    // Private data

        //- Name of volume fraction field to measure
        word vofField_;

        scalar contourLevel_;

        //- Names of gauges from input dict
        wordList gaugeNames_;

        //- Output file header information
        virtual void writeFileHeader(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("waveGauge");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        //  allow the possibility to load fields from files
        waveGauge
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        waveGauge(const waveGauge&) = delete;


    //- Destructor
    virtual ~waveGauge();


    // Member Functions

        //- Read the waveGauge
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const waveGauge&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
