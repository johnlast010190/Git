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
    (c) 2024-2025 Engys Ltd.

Class
    Foam::fvMeshGIBChangers::patchSwitch

Description
    Dynamic GIB mesh changer for patch switches.

SourceFiles
    fvMeshGIBChangersPatchSwitch.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshGIBChangersPatchSwitch_H
#define fvMeshGIBChangersPatchSwitch_H

#include "fvMesh/fvMeshGIBChangers/fvMeshGIBChanger/fvMeshGIBChanger.H"
#include "patchSwitch/switches/GIBSwitch/GIBSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{

/*---------------------------------------------------------------------------*\
                Class fvMeshGIBChangers::patchSwitch Declaration
\*---------------------------------------------------------------------------*/

class patchSwitch
:
    public fvMeshGIBChanger
{
    // Private Data

        //- Dynamic mesh dictionary
        dictionary dynamicMeshCoeffs_;

        //- Switches to activate patches
        PtrList<GIBSwitch> pSwitchPtr_;


public:

    //- Runtime type information
    TypeName("patchSwitch");


    // Constructors

        //- Construct from IOobject
        explicit patchSwitch(fvMesh&);

        //- Construct from fvMesh and dictionary
        patchSwitch(fvMesh&, const dictionary);

        //- Disallow default bitwise copy construction
        patchSwitch(const patchSwitch&) = delete;


    //- Destructor
    ~patchSwitch(){};


    // Member Functions

        //- Update the mesh for changes in GIB patches
        virtual bool update();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const patchSwitch&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshGIBChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
