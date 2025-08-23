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
    TimelineCollection

Description
    Collects a HashTable of Timelines.

    If no timeline is currently present then one is added with a header according
    to the specification

SourceFiles
    TimelineCollection.C

Contributors/Copyright:
    2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef TimelineCollection_H
#define TimelineCollection_H

#include "db/objectRegistry/objectRegistry.H"
#include "include/swakTime.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class TimelineCollection Declaration
\*---------------------------------------------------------------------------*/

class TimelineCollection
{
    // Private data

        //- The output-directory
        fileName outputDirectory_;

        //- The time object
        const Time& time_;

        //- Header specifications
        dictionary headerSpecs_;

        //- Output file pointers
        HashPtrTable<OFstream> outputFilePtr_;

        //- Default time name
        word timeName_;


public:

    // Static data members

    TypeName("TimelineCollection");

    // Constructors

        //- Construct null
        TimelineCollection
        (
            const fileName& dir,
            const Time& time,
            bool useStartTime = true
        );

        //- Disallow default bitwise copy construction
        TimelineCollection(const TimelineCollection&) = delete;


    //- Destructor
    virtual ~TimelineCollection();


    // Member Functions

        void addSpec
        (
            const string& nameExpr,
            const string& headerSpecString
        );


    // Member Operators

        OFstream& operator()
        (
            const word& name,
            bool autoParallel = true
        );

        OFstream& operator[](const word& name)
        {
            return (*this)(name, false);
        }

        //- Disallow default bitwise assignment
        void operator=(const TimelineCollection&) = delete;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
