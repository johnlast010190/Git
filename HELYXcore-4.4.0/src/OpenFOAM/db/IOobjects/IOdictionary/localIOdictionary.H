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
    (c) 2015-2020 OpenFOAM Foundation

Class
    Foam::localIOdictionary

Description
    localIOdictionary is derived from IOdictionary with global flag set
    to false to disable parallel master reading.

    Used to read GeometricFields.

SourceFiles
    localIOdictionary.C

\*---------------------------------------------------------------------------*/

#ifndef localIOdictionary_H
#define localIOdictionary_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class localIOdictionary Declaration
\*---------------------------------------------------------------------------*/

class localIOdictionary
:
    public IOdictionary
{
public:

    // Constructors

        //- Construct given an IOobject
        localIOdictionary(const IOobject& io);

        //- Construct given an IOobject and dictionary
        localIOdictionary(const IOobject&, const dictionary&);

        //- Construct given an IOobject and actual type name
        localIOdictionary(const IOobject& io, const word& actualType);


    // Member functions

        //- Is object global
        virtual bool global() const
        {
            return false;
        }

        //- Return complete path + object name if the file exists
        //  in the case otherwise null
        virtual fileName filePath() const
        {
            // Use default (local only) search strategy
            return localFilePath(type());
        }
};


//- Template function for obtaining global status
template<>
inline bool typeGlobal<localIOdictionary>()
{
    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
