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
    surfaceOptimizer

Description
    Performs optimisation of a central vertex in a simplex

SourceFiles
    surfaceOptimizer.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceOptimizer_H
#define surfaceOptimizer_H

#include "meshes/primitiveShapes/point/point.H"
#include "meshes/meshShapes/triFace/triFace.H"
#include "utilities/containers/DynList/DynList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class surfaceOptimizer Declaration
\*---------------------------------------------------------------------------*/

class surfaceOptimizer
{
    // Private static data

        //- Direction vectors for divide and conquer algorithm
        static const vector dirVecs[4];


    // Private data

        //- Reference to the simplex points
        DynList<point>& pts_;

        //- Reference to the triangles forming a simplex
        const DynList<triFace>& trias_;

        //- Min position of the bnd box
        point pMin_;

        //- Max position of the bnd box
        point pMax_;


    // Private member functions

        //- Evaluate stabilisation factor
        scalar evaluateStabilisationFactor() const;

        //- Evaluate the functional
        scalar evaluateFunc(const scalar& K) const;

        //- Evaluate gradients needed for optimisation
        void evaluateGradients(const scalar&, vector&, tensor&) const;

        //- Optimise point position using the divide and conquer technique
        scalar optimiseDivideAndConquer(const scalar tol);

        //- Optimise point position via the steepest descent method
        scalar optimiseSteepestDescent(const scalar tol);


public:

    // Constructors

        //- Construct from transformed points and triangles forming a simplex
        surfaceOptimizer
        (
            DynList<point>& pts,
            const DynList<triFace>& trias
        );

        //- Disallow default bitwise copy construction
        surfaceOptimizer(const surfaceOptimizer&) = delete;


    //- Destructor
    ~surfaceOptimizer();


    // Member Functions

        //- Optimizes position of a central point in the simplex
        point optimizePoint(const scalar tol = 0.1);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const surfaceOptimizer&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
