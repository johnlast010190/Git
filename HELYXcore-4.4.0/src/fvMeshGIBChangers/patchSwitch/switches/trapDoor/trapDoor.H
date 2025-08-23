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
    (c) 2019-2024 Engys Ltd.

Class
    Foam::trapDoor

Description
    Active patch based on pressure drop. Deactivate patch based on normal
    velocity of fluid.

    If the alpha field does not exist, take the velocity of the whole patch
    into accoun,t which will work for single phase flow.

SourceFiles
    trapDoor.C

\*---------------------------------------------------------------------------*/

#ifndef trapDoor_H
#define trapDoor_H

#include "patchSwitch/switches/GIBSwitch/GIBSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{

/*---------------------------------------------------------------------------*\
                           Class trapDoor Declaration
\*---------------------------------------------------------------------------*/

class trapDoor
:
    public GIBSwitch
{
    // Private Member Functions

        //- Direction invertion
        Switch invertDirection_;

        //- Velocity threshold
        scalar fVel_;

        //- Pressure drop threshold
        scalar dp_;

        //- Name of velocity field
        word UName_;

        //- Name of pressure field
        word pName_;

        //- Name of fraction field
        word fractionName_;

        virtual bool enableCondition();

        virtual bool disableCondition();


public:

    //- Runtime type information
    TypeName("trapDoor");


    // Constructors

        //- Construct from fvMesh and dictionary
        trapDoor(const fvMesh&, const dictionary&);

        //- Disallow default bitwise copy construction
        trapDoor(const trapDoor&) = delete;


    //- Destructor
    virtual ~trapDoor(){};


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const trapDoor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshGIBChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
