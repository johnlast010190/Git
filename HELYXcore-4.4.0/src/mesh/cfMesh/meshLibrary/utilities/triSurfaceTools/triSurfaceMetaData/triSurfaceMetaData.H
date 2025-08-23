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
    triSurfaceMetaData

Description
    Generates meta data of the surface mesh. It contains:
    The number of points
    The number of triangles
    The number of patches
    The number of feature edges
    The number of subsets, names, and the number of elements in each subset

SourceFiles
    triSurfaceMetaData.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceMetaData_H
#define triSurfaceMetaData_H

#include "utilities/meshes/triSurf/triSurf.H"
#include "db/dictionary/dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class triSurfaceMetaData Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceMetaData
{
    // Private data

        //- Const reference to triSurf
        const triSurf& surf_;

        //- Patches/subsets for removal
        dictionary metaDict_;


    // Private member functions

        //- Create dictionary with meta data
        void createMetaData();


public:

    // Constructors

        //- Construct from triSurf
        triSurfaceMetaData(const triSurf& surface);

        //- Disallow default bitwise copy construction
        triSurfaceMetaData(const triSurfaceMetaData&) = delete;


    //- Destructor
    ~triSurfaceMetaData();


    // Member Functions

        //- Return a constant reference to meta data
        const dictionary& metaData() const
        {
            return metaDict_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceMetaData&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
