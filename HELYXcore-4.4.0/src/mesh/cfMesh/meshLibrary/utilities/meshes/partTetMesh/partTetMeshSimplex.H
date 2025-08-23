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

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Class
    partTetMeshSimplex

Description
    A simplex which is used for smoothing purposes

SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef partTetMeshSimplex_H
#define partTetMeshSimplex_H

#include "utilities/meshes/partTetMesh/partTetMesh.H"
#include "utilities/helperClasses/parallelHelpers/parPartTet/parPartTet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class partTetMesh;

/*---------------------------------------------------------------------------*\
                           Class partTetMeshSimplex Declaration
\*---------------------------------------------------------------------------*/

class partTetMeshSimplex
{
    // Private data
        //- points making the simplex
        DynList<point, 128> pts_;

        //- tets making the simplex
        DynList<partTet, 128> tets_;

    public:

    // Constructors
        //- Construct from partTetMeshSimplex and point label
        partTetMeshSimplex(const partTetMesh& tm, const label pI);

        //- Construct from the list of parPartTet and point label
        partTetMeshSimplex
        (
            const DynList<parPartTet>& pt,
            const label gpI
        );

        //- contruct from a list of points and a list of partTets
        //- ad a point label
        partTetMeshSimplex
        (
            const DynList<point, 128>& pts,
            const DynList<partTet, 128>& tets,
            const label pointI
        );

    // Destructor
        ~partTetMeshSimplex();

    // Member functions
        //- return points
        inline DynList<point, 128>& pts()
        {
            return pts_;
        }

        //- return points
        inline const DynList<point, 128>& pts() const
        {
            return pts_;
        }

        //- return tets
        inline const DynList<partTet, 128>& tets() const
        {
            return tets_;
        }

        //- return centre point coordinates
        inline const point& centrePoint() const
        {
            return pts_[tets_[0][3]];
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
