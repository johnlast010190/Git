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
    TimeClone

Description
    Object that on request creates a full copy of all fields etc and writes them on
    request

    Only copies objectRegistries or objects that are auto-writeable

SourceFiles
    TimeClone.C

Contributors/Copyright:
    2014-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef TimeClone_H
#define TimeClone_H

#include "db/objectRegistry/objectRegistry.H"
#include "include/swakTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class TimeClone Declaration
\*---------------------------------------------------------------------------*/

class TimeClone
{
    // Private data

        //- The actual time-thing
        autoPtr<Time> storedTime_;


    // Private Member Functions

        //- Do the actual copying
        label copyObjects(const objectRegistry& src, objectRegistry& dst);


public:

    //- Runtime type information
    TypeName("TimeClone");

    // Constructors

        //- Construct null
        TimeClone();

        //- Disallow default bitwise copy construction
        TimeClone(const TimeClone&) = delete;


    // Selectors

        //- Make a copy
        void copy(const Time&);


    //- Destructor
    virtual ~TimeClone();


    // Member Functions

        //- Write the time-directory
        bool write(const bool force = false);

        bool ok() const
        {
            return storedTime_.valid();
        }


    // Member Operators

        const Time& operator()() const;

        //- Disallow default bitwise assignment
        void operator=(const TimeClone&) = delete;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
