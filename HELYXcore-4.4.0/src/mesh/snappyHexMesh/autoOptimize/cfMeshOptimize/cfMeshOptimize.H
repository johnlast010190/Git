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
    (c) 2014-1014 Engys Ltd.

Class
    Foam::cfMeshOptimize

Description
    Perform optimisation of mesh

SourceFiles
    cfMeshOptimize.C

\*---------------------------------------------------------------------------*/

#ifndef cfMeshOptimize_H
#define cfMeshOptimize_H

#include "fvMesh/fvMesh.H"
#include "autoOptimize/autoOptimize.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class cfMeshOptimize Declaration
\*---------------------------------------------------------------------------*/

class cfMeshOptimize : public autoOptimize
{

private:

      // Private data

        //- Mesh
        fvMesh& mesh_;

        //- Relative threshold for worst quality faces
        scalar qualityThreshold_;

        //- Number of Laplacian smoothing iterations
        label numLaplaceIterations_;

        //- Maximum global number of optimisation iteration
        label maxNumGlobalIterations_;

        //- Maximum number local optimisation iteration
        label maxNumIterations_;

        //- Maximum number of surface optimisation iteration
        label maxNumSurfaceIterations_;

        //- If set > -1 only optimise error cells and specified buffer layers
        label errorBufferLayers_;

        //- Whether to perform a set of relaxed quality checks
        Switch relaxedCheck_;

        //- Whether to perform a set of relaxed quality checks for
        // boundary motion
        Switch relaxedBoundaryCheck_;

        //- Whether to perform checks on warped faces
        Switch checkWarped_;

        //- Whether to perform final sync of optimised points
        Switch syncPts_;

        //- Whether to perform orthogonality optimisation (180 disable)
        scalar maxOrth_;

        //- Tolerance used for minimum face areas
        scalar minFaceArea_;

public:

    //- Runtime type information
    TypeName("cfMeshOptimize");


        //- Update input points field
        void movePoints(pointField& newPoints);

        //- optimise mesh points
        void optimize();


        // Constructors
        cfMeshOptimize
        (
            fvMesh& mesh,
            const dictionary& dict
        );

        //- Destructor
        virtual ~cfMeshOptimize();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
