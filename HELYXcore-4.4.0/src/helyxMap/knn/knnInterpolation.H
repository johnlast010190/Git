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
    (c) 2024-2025 Engys Ltd.

Class
    Foam::knnInterpolation

Description:
    K-nearest neighbour interpolation class.

SourceFiles
    knnInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef knnInterpolation_H
#define knnInterpolation_H

#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class knnInterpolation Declaration
\*---------------------------------------------------------------------------*/

class knnInterpolation
{
protected:

    // Protected Data

        const label numNbrs_;

        bool useInvDistSq_;

        autoPtr<pointField> sourcePoints_;


    // Protected Member Functions

        //- Interpolate from k nearest neighbours of sourceV vector
        //  using inverse distance
        point interpolateVector
        (
            const std::vector<std::vector<scalar>>& sourceV,
            const std::vector<label>& nb,
            const std::vector<scalar>& dist
        );

        //- Interpolate from k nearest neighbours of sourceV,
        //  but weight surface area as well.
        point interpolateVector
        (
            const std::vector<std::vector<scalar>>& sourceV,
            const scalarField& sourceAreas,
            const std::vector<label>& nb,
            const std::vector<scalar>& dist
        );

        //- Interpolate from k nearest neighbours of sourceV,
        //  but weight with orthogonality as well.
        point interpolateVector
        (
            const std::vector<std::vector<scalar>>& sourceV,
            const pointField& sourceNormal,
            const point& ptNormal,
            const std::vector<label>& nb,
            const std::vector<scalar>& dist
        );

        //- Interpolate from k nearest neighbours of sourceV
        //  using inverse distance
        scalar interpolateScalar
        (
            const std::vector<scalar>& sourceV,
            const std::vector<label>& nb,
            const std::vector<scalar>& dist
        );

        //- For each point of the pointField find the closest triangle (1)
        // from triSurfaceMesh; nearest function of triSurfaceMesh is used.
        scalar calculateDistance
        (
            const pointField& actualPnt,
            const triSurfaceMesh& targetSM
        );


public:

    // Static Member Functions

        static wordList timeNames(const instantList& times);

        static bool findTime
        (
            const instantList& times,
            const label startSampleTime,
            const scalar timeVal,
            label& lo,
            label& hi
        );


    // Constructors

        knnInterpolation
        (
            const label numNeighbours,
            bool useInverseDistanceSquared,
            const pointField& inputPoints
        );

        knnInterpolation
        (
            const label numNeighbours,
            bool useInverseDistanceSquared,
            const fileName& pointsFile
        );

        knnInterpolation
        (
            const label numNeighbours,
            bool useInverseDistanceSquared
        );


    // Member Functions

        template <class Type>
        void field2Vector
        (
            const Type& vf,
            std::vector<std::vector<scalar>>& vect
        ) const
        {
            forAll(vf, pi)
            {
                std::vector<scalar> pt(3);
                pt[0] = vf[pi].x();
                pt[1] = vf[pi].y();
                pt[2] = vf[pi].z();

                vect.push_back(pt);
            }
        }

        template <class Type>
        void field2Scalar(const Type& vf, std::vector<scalar>& vect) const
        {
            forAll(vf, pi)
            {
                vect.push_back(vf[pi]);
            }
        }

        tmp<scalarField> calculateInterpolatedScalar
        (
            const scalarField& sourceField,
            const vectorField& faceCentres
        );

        tmp<vectorField> calculateInterpolatedVector
        (
            const vectorField& sourceField,
            const vectorField& faceCentres
        );

        tmp<pointField> readFile(const fileName& fn);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
