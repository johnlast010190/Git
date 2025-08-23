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
    (c) 2024 Engys Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Class
    meshSurfaceEngineModifier

Description
    Modifier for the meshSurfaceEngine

SourceFiles
    meshSurfaceEngine.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceEngineModifier_H
#define meshSurfaceEngineModifier_H

#include "utilities/surfaceTools/meshSurfaceEngine/meshSurfaceEngine.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class meshSurfaceEngineModifier Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceEngineModifier
{
    // Private data

        //- Reference to the meshSurfaceEngine
        meshSurfaceEngine& surfaceEngine_;


public:

    // Constructors

        //- Construct from reference to meshSurfaceEngine
        meshSurfaceEngineModifier(meshSurfaceEngine& surfaceEngine);

        //- Construct from const reference to meshSurfaceEngine
        meshSurfaceEngineModifier(const meshSurfaceEngine& surfaceEngine);

        //- Disallow default bitwise copy construction
        meshSurfaceEngineModifier(const meshSurfaceEngineModifier&) = delete;


    //- Destructor
    ~meshSurfaceEngineModifier();


    // Member Functions

        //- Relocate the selected boundary vertex
        void moveBoundaryVertexNoUpdate(const label bpI, const point& newP);

        //- Relocate the selected boundary vertex and update geometry data
        void moveBoundaryVertex(const label bpI, const point& newP);

        //- Update normals of boundary vertices at processor boundaries
        inline void updateVertexNormals()
        {
            if (surfaceEngine_.pointNormalsPtr_)
            {
                surfaceEngine_.updatePointNormalsAtProcBoundaries();
            }
        }

        //- Updates faceCentres, face normal, and point normals
        //  of the boundary faces attached to the selected boundary points
        void updateGeometry(const labelLongList&);
        void updateGeometry();

        //- Makes sure that all surface vertices at parallel boundaries
        //  have the same coordinates
        void syncVerticesAtParallelBoundaries();
        void syncVerticesAtParallelBoundaries(const labelLongList&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceEngineModifier&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
