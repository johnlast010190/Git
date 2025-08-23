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
    Foam::pointProfileDataList

Description
    Container for multiple pointProfileData sets

SourceFiles
    pointProfileDataList.C

\*---------------------------------------------------------------------------*/

#ifndef pointProfileDataList_H
#define pointProfileDataList_H

#include "fields/fvPatchFields/derived/mappedProfile/pointProfileData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class pointProfileDataList Declaration
\*---------------------------------------------------------------------------*/

class pointProfileDataList
:
    public PtrList<pointProfileData>
{
public:

    // Constructors

        //- Construct null
        pointProfileDataList();

        //- Construct from Istream
        pointProfileDataList(ISstream&);

        //- Disallow default bitwise copy construction
        pointProfileDataList(const pointProfileDataList&) = delete;


    //- Destructor
    ~pointProfileDataList();


    // Member Functions

        //- Provide mapped data
        scalarField* mapData(const pointField&, word, word);

        //- Read
        void read(ISstream&);


    // Member Operators


        //- Disallow default bitwise assignment
        void operator=(const pointProfileDataList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
