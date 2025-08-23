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
    (c) 2014-2022 OpenFOAM Foundation

Class
    Foam::IOdistributionMap

Description
    IOdistributionMap is derived from distributionMap and
    IOobject to give the distributionMap
    automatic IO functionality via the objectRegistry.

SourceFiles
    IOdistributionMap.C

\*---------------------------------------------------------------------------*/

#ifndef IOdistributionMap_H
#define IOdistributionMap_H

#include "meshes/polyMesh/polyDistributionMap/distributionMap.H"
#include "db/regIOobject/regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class IOdistributionMap Declaration
\*---------------------------------------------------------------------------*/

class IOdistributionMap
:
    public regIOobject,
    public distributionMap
{

public:

    //- Runtime type information
    TypeName("distributionMap");

    // Constructors

        //- Construct given an IOobject
        IOdistributionMap(const IOobject&);

        //- Construct given an IOobject and distributionMap
        IOdistributionMap(const IOobject&, const distributionMap&);

        //- Move constructor by transferring the distributionMap contents
        IOdistributionMap(const IOobject&, distributionMap&&);


    //- Destructor
    virtual ~IOdistributionMap();


    // Member functions

        //- ReadData function required for regIOobject read operation
        virtual bool readData(Istream&);

        //- WriteData function required for regIOobject write operation
        virtual bool writeData(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
