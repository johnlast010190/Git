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
    triSurfaceCurvatureCalculator

Description
    Provides information regarding surface partitions on the surface
    triangulation, and connectivity between various surface partitions.

SourceFiles
    triSurfaceCurvatureEstimator.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceCurvatureEstimator_H
#define triSurfaceCurvatureEstimator_H

#include "utilities/meshes/triSurf/triSurf.H"
#include "utilities/containers/FRWGraph/FRWGraph.H"
#include "utilities/containers/DynList/DynList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class triSurfaceCurvatureEstimator Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceCurvatureEstimator
{
    // Private data

        //- Reference to triSurface
        const triSurf& surface_;

        //- Curvature of points at feature edges
        scalarField edgePointCurvature_;

        //- Curvatures of other points
        FRWGraph<label, 3> patchPositions_;
        List<DynList<scalar, 1>> gaussianCurvature_;
        List<DynList<scalar, 1>> meanCurvature_;
        List<DynList<scalar, 1>> maxCurvature_;
        List<DynList<scalar, 1>> minCurvature_;
        List<DynList<vector, 1>> maxCurvatureVector_;
        List<DynList<vector, 1>> minCurvatureVector_;


    // Private member functions

        //- Calculate curvature of feature edges
        void calculateEdgeCurvature();

        //- Calculate curvatures of other surface points
        void calculateSurfaceCurvatures();

        //- Calculate Gaussian curvature
        void calculateGaussianCurvature();

        //- Calculation of mean curvature
        void calculateMeanCurvature();

        //- Calculation of min and max curvature
        void calculateMinAndMaxCurvature();

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceCurvatureEstimator&) = delete;


public:

    // Constructors

        //- Construct from triSurface
        triSurfaceCurvatureEstimator(const triSurf& surface);

        //- Disallow default bitwise copy construction
        triSurfaceCurvatureEstimator
        (
            const triSurfaceCurvatureEstimator&
        ) = delete;


    //- Destructor
    ~triSurfaceCurvatureEstimator();


    // Member Functions

        //- Return curvature of points at feature edges
        scalar edgePointCurvature(const label) const;
        scalar curvatureAtEdge(const label) const;

        //- Return average curvature on a surface triangle
        scalar gaussianCurvatureAtTriangle(const label) const;
        scalar meanCurvatureAtTriangle(const label) const;
        scalar maxCurvatureAtTriangle(const label) const;
        scalar minCurvatureAtTriangle(const label) const;

        vector maxCurvatureVectorAtTriangle(const label) const;
        vector minCurvatureVectorAtTriangle(const label) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
