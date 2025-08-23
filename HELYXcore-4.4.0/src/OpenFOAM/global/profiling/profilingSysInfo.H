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
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::profilingSysInfo

Description
    General system information useful for profiling

SourceFiles
    profilingSysInfo.C

\*---------------------------------------------------------------------------*/

#ifndef profilingSysInfo_H
#define profilingSysInfo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Ostream;
class profilingSysInfo;

/*---------------------------------------------------------------------------*\
                      Class profilingSysInfo Declaration
\*---------------------------------------------------------------------------*/

class profilingSysInfo
{
public:


    // Constructors

        //- Construct from components
        profilingSysInfo();

        //- Disallow default bitwise copy construction
        profilingSysInfo(const profilingSysInfo&) = delete;


    //- Destructor
    ~profilingSysInfo();


    // Member Functions

        //- Update it with a new timing information
        void update();

        //- Write the profiling system-info, use dictionary format.
        Ostream& write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const profilingSysInfo&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
