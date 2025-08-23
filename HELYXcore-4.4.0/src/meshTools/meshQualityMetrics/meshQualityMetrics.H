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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::meshQualityMetrics

Description
    MeshObject that marks poor quality cells (faces) that can be used from
    other classes. Currently four checks are performed:
    -cell volumes
    -face areas
    -face orthogonality
    -face pyramids

SourceFiles
    meshQualityMetrics.C

\*---------------------------------------------------------------------------*/

#ifndef meshQualityMetrics_H
#define meshQualityMetrics_H

#include "meshes/MeshObject/MeshObject.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "primitives/ints/lists/labelList.H"
#include "primitives/Vector/vector/vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                     Class meshQualityMetrics Declaration
\*---------------------------------------------------------------------------*/

class meshQualityMetrics
:
    public MeshObject<polyMesh, UpdateableMeshObject, meshQualityMetrics>
{
    // Private data

        //- Const dictionary

        //- Marker for bad quality cells based on criteria
        mutable boolList* badQualityCellsPtr_;

        //- MeshObject dictionary
        dictionary mqDict_;

        //- Report statistics from mesh quality functions
        Switch reportStatistics_;

        //- Write fields on disc for debugging
        Switch writeField_;

        //- Min volume threshold for face pyramid check
        scalar minVol_;

        //- Min area threshold
        scalar minArea_;

        //- Maximum allowed non-orthogonality threshold
        scalar maxNonOrtho_;

        //- Maximum allowed skewness threshold
        scalar maxSkew_;

        //- Min determinant
        scalar minDeterminant_;

        //- Min faceWeight
        scalar minFaceWeight_;

    // Private Member Functions

        //- read thresholds
        void readMetricThresholds();

        //- Mark bad cells
        void calcMetrics() const;

        //- Clear list
        void clearMetrics() const;

public:

    // Declare name of the class and its debug switch
    ClassName("meshQualityMetrics");


    // Constructors

        //- Construct from components
        meshQualityMetrics(const polyMesh&, const dictionary&, const word&);

        //- Disallow default bitwise copy construction
        meshQualityMetrics(const meshQualityMetrics&) = delete;


    //- Destructor
    ~meshQualityMetrics();

    // Selectors

        //- Constructor to create new meshQualityMetric based on meshObject
        // dictionary in system
        static meshQualityMetrics& New
        (
            const polyMesh&,
            const word&
        );

    // Member Functions

        //- Update topology
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);

        //- Correct weighting factors for moving mesh.
        virtual bool movePoints();

        // Access

        //- Access the bad cell marker
        const boolList& badQualityCells() const;

        //- write field to disc for debudding
        const Switch& writeField() const
        {
            return writeField_;
        }

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshQualityMetrics&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
