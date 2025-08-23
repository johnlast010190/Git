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
    triSurface2DCheck

Description


SourceFiles
    triSurface2DCheck.C

\*---------------------------------------------------------------------------*/

#ifndef triSurface2DCheck_H
#define triSurface2DCheck_H

#include "utilities/meshes/triSurf/triSurf.H"
#include "primitives/SymmTensor/symmTensor/symmTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class triSurface2DCheck Declaration
\*---------------------------------------------------------------------------*/

class triSurface2DCheck
{
    // Private data

        //- Reference to triSurf
        const triSurf& surf_;

        //- Covariance matrix
        symmTensor covarianceMatrix_;


    // Private member functions

        //- Create covariance matrix
        void createCovarianceMatrix();


public:

    // Constructors

        //- Construct from octree
        triSurface2DCheck(const triSurf& surface);

        //- Disallow default bitwise copy construction
        triSurface2DCheck(const triSurface2DCheck&) = delete;


    //- Destructor
    ~triSurface2DCheck();


    // Member Functions

        //- checks if the surface is a 2D triangulation
        bool is2DSurface() const;

        //- create subset containing invalid facets
        void createSubsets();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const triSurface2DCheck&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
