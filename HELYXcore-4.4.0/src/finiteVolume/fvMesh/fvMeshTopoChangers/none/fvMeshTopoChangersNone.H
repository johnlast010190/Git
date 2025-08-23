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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2023 Engys Ltd.

Class
    Foam::fvMeshTopoChangers::none

Description
    Dummy fvMeshTopoChanger which does not apply any topology change.

SourceFiles
    none.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshTopoChangersNone_H
#define fvMeshTopoChangersNone_H

#include "fvMesh/fvMeshTopoChangers/fvMeshTopoChanger/fvMeshTopoChanger.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{

/*---------------------------------------------------------------------------*\
                   Class fvMeshTopoChangers::none Declaration
\*---------------------------------------------------------------------------*/

class none
:
    public fvMeshTopoChanger
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from fvMesh
        none(fvMesh& mesh);

        //- Construct from fvMesh and dictionary
        none(fvMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        none(const none&) = delete;


    //- Destructor
    ~none();


    // Member Functions

        //- The mesh is static and does not change, so return false
        virtual bool dynamic() const override
        {
            return false;
        }

        //- Dummy update function which does not change the mesh
        virtual bool update() override
        {
            return false;
        }

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&) override
        {}

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&) override
        {}

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&) override
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const none&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshTopoChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
