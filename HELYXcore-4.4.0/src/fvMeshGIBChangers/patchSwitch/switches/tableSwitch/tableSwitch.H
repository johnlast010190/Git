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
    Foam::tableSwitch

Description
    Activate/Deactivate based on user-defined table

SourceFiles
    tableSwitch.C

\*---------------------------------------------------------------------------*/

#ifndef tableSwitch_H
#define tableSwitch_H

#include "patchSwitch/switches/GIBSwitch/GIBSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{

/*---------------------------------------------------------------------------*\
                         Class tableSwitch Declaration
\*---------------------------------------------------------------------------*/

class tableSwitch
:
    public GIBSwitch
{
    // Private Data

        //- Switch table
        autoPtr<Function1<bool>> boolTable_;


    // Private Member Functions

        virtual bool enableCondition();

        virtual bool disableCondition();


public:

    //- Runtime type information
    TypeName("tableSwitch");


    // Constructors

        //- Construct from fvMesh and dictionary
        tableSwitch(const fvMesh&, const dictionary&);

        //- Disallow default bitwise copy construction
        tableSwitch(const tableSwitch&) = delete;


    //- Destructor
    virtual ~tableSwitch(){};


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const tableSwitch&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshGIBChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
