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
    TimeCloneList

Description
    List of N time-clones. The oldest is always removed. After writing all are removed

    Uses a list of Ptr because PtrList does not allow the direct manipulation of the pointer values
SourceFiles
    TimeCloneList.C

Contributors/Copyright:
    2014-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef TimeCloneList_H
#define TimeCloneList_H

#include "TimeClone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class TimeCloneList Declaration
\*---------------------------------------------------------------------------*/

class TimeCloneList
{
    // Private data

        //- the actual time-thing
        List<TimeClone*> storedTimes_;

        static label count_;


    // Private Member Functions

        void clear();


public:

    //- Runtime type information
    TypeName("TimeCloneList");


    // Constructors

        //- Construct null
        TimeCloneList(const dictionary& dict);

        //- Disallow default bitwise copy construction
        TimeCloneList(const TimeCloneList&) = delete;


    // Selectors

        //- Make a copy
        void copy(const Time&);


    //- Destructor
    virtual ~TimeCloneList();


    // Member Functions

        //- Write the time-directory
        bool write(const bool force = false);

        bool has(const Time& other);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const TimeCloneList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
