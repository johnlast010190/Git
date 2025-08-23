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
    (c) 2022 OpenFOAM Foundation
    (c) 2024-2025 Engys Ltd.

Class
    Foam::fvMeshTopoChangers::meshToMesh

Description
    fvMeshTopoChanger which maps to new mesh

SourceFiles
    fvMeshTopoChangersMeshToMesh.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshTopoChangersMeshToMesh_H
#define fvMeshTopoChangersMeshToMesh_H

#include "fvMesh/fvMeshTopoChangers/fvMeshTopoChanger/fvMeshTopoChanger.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{

/*---------------------------------------------------------------------------*\
                Class fvMeshTopoChangers::meshToMesh Declaration
\*---------------------------------------------------------------------------*/

class meshToMesh
:
    public fvMeshTopoChanger
{
    // Private Data

        //- Mesh-to-mesh control dictionary
        dictionary dict_;

        //- List of mesh mapping times
        scalarList times_;

        //- Time delta used for time -> index
        scalar timeDelta_;

        //- Hash set of mesh mapping time indices
        labelHashSet timeIndices_;

        //- The time index used for updating
        label timeIndex_;


    // Private Member Functions

        //- Interpolate U's to Uf's
        void interpolateUfs();


public:

    //- Runtime type information
    TypeName("meshToMesh");


    // Constructors

        //- Construct from fvMesh and dictionary
        meshToMesh(fvMesh&, const dictionary& dict);

        //- Disallow default bitwise copy construction
        meshToMesh(const meshToMesh&) = delete;


    //- Destructor
    virtual ~meshToMesh();


    // Member Functions

        const scalarList& times() const
        {
            return times_;
        }

        scalar timeDelta() const
        {
            return timeDelta_;
        }

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshToMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshTopoChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
