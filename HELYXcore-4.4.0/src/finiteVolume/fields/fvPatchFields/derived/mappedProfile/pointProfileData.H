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
    (c) 2009-2009 OpenCFD Ltd.
    (c) 2017-2024 Engys Ltd.

Class
    Foam::pointProfileData

Description
    Class for readfing and mapping fluent profile input files to patches
    Currently only for scalars

SourceFiles
    pointProfileData.C

\*---------------------------------------------------------------------------*/

#ifndef pointProfileData_H
#define pointProfileData_H

#include "db/dictionary/dictionary.H"
#include "db/IOstreams/IOstreams/Istream.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "fields/Fields/scalarField/scalarField.H"
#include "meshes/primitiveShapes/point/pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class pointProfileData Declaration
\*---------------------------------------------------------------------------*/

class pointProfileData
{
    // Private data

        //- boundary region name
        word regionName_;

        //- size
        label N_;

        //- field names
        wordList fieldNames_;

        //- fields
        PtrList<scalarField> fields_;

        //- point positions (derived from xyz fields)
        pointField points_;


public:

    // Constructors

        //- Construct null
        pointProfileData();

        //- Construct from Istream
        pointProfileData(ISstream&);

        //- Disallow default bitwise copy construction
        pointProfileData(const pointProfileData&) = delete;


    //- Destructor
    ~pointProfileData();


    // Member Functions

        //- Name of the profile boundary region
        const string& regionName() const
        {
            return regionName_;
        }

        //- Map requested field to provided points
        scalarField* mapData(const pointField&, word);

        //- Read
        void read(ISstream&);

        // unexpected character read error
        static void checkIO(const token&, const token::punctuationToken&);

        // go to next valid token
        static void nextValidToken(ISstream&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pointProfileData&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
