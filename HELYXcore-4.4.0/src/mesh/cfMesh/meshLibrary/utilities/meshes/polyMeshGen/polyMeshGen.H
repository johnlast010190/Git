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
    (c) Creative Fields, Ltd.
    (c) 2017 Engys Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Class
    polyMeshGen

Description
    Mesh with selections

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGen_H
#define polyMeshGen_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "utilities/meshes/polyMeshGen/polyMeshGenCells.H"

#include "db/dictionary/dictionary.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMeshGen
: public polyMeshGenCells
{
    // Private data
        //- meta data about the meshing process
        IOdictionary metaDict_;

public:

    // Constructors
        //- Null constructor
        polyMeshGen(const Time& t);

        //- Construct from components without the boundary
        polyMeshGen
        (
            const Time& t,
            const pointField& points,
            const faceList& faces,
            const cellList& cells
        );

        //- Construct from components with the boundary
        polyMeshGen
        (
            const Time& t,
            const pointField& points,
            const faceList& faces,
            const cellList& cells,
            const polyBoundaryMesh& patches
        );

    // Destructor

        ~polyMeshGen();

    // Public member functions
        //- return a constant reference to metaDict
        inline const dictionary& metaData() const
        {
            return metaDict_;
        }

        //- return a reference to metaDict
        inline dictionary& metaData()
        {
            return metaDict_;
        }

    // Read mesh
        void read();

    // Write mesh
        void write() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
