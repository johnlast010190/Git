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
    (c) 2020-2024 Engys Ltd

Class
    Foam::VDBGridsConverter

Description
    Converts a std:vector of VDB grids into an fvMesh

\*---------------------------------------------------------------------------*/

#ifndef VDBGridsConverter_H
#define VDBGridsConverter_H

#include "helyxVDB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class VDBGridsConverter Declaration
\*---------------------------------------------------------------------------*/


class VDBGridsConverter
{
    using gridLists = std::vector<std::vector<FloatGrid::Ptr>>;

    // Private data

        //- Reference to the helyxVDB object
        helyxVDB& hVDB_;

        //- Reference to the VDB grids to be converted
        gridLists& allCellsGrids_;

        //-
        label pointsPerEdge_; //TODO -> edgeDivisions

        //- Point list in index space
        pointList xyzPointsList; //TODO -> xyzPoints_

        //- Boundary info
        List<word>  boundaryPatchNames_;
        List<word>  boundaryPatchTypes_;
        List<label> boundaryPatchStarts_;
        List<label> boundaryPatchSizes_;

        //- Result points
        pointField points_;

        //- Result faces
        faceList faces_;

        //- Result owner
        labelList owner_;

        //- Result neighbour
        labelList neighbour_;

        //- Result fields
        scalarField cellLevelField_;
        scalarField distanceField_;
        scalarField patchIDField_;
        scalarField procIDField_;


    // Private Member Functions

        //- Returns coordinates of hex associated with voxel
        static const std::vector<openvdb::Coord> getVoxelVertices
        (
            const openvdb::Coord& coord,
            const label pointsPerEdge = 2
        );

        //- Returns coordinates of hex vertices associated with voxel
        const std::vector<openvdb::Coord> getVoxelFaceCentres
        (
            const openvdb::Coord& coord,
            const label pointsPerEdge = 2
        );

        //- true if a coarse face at coord xyx overlaps with 4 faces
        // of finer level
        struct InterfaceFace
        {
            bool isInterfaceFace;
            openvdb::Coord xyz;
        };

        //- Finds out if a face overlaps with 4 of finer level
        // (therefore the cell needs to change from hex to slitHex)
        template<typename GridT>
        InterfaceFace isCellLevelInterfaceFace
        (
            const face& face,
            openvdb::math::Transform::ConstPtr coarseTransform,
            typename GridT::ConstAccessor& cCoarseFaceGridAccessor
        );

        //-
        void convert();


public:

    // Constructors

        //- Construct from helyxVDB and lists of VDBGrids
        VDBGridsConverter
        (
            helyxVDB&,
            gridLists&,
            label pointsPerEdge = 2
        );

        //- Disallow default bitwise copy construction
        VDBGridsConverter(const VDBGridsConverter&) = delete;


    //- Destructor
    ~VDBGridsConverter();


    // Member Functions

        // Access

        List<word>&  boundaryPatchNames()
        {
            return boundaryPatchNames_;
        }

        List<word>&  boundaryPatchTypes()
        {
            return boundaryPatchTypes_;
        }

        List<label>& boundaryPatchStarts()
        {
            return boundaryPatchStarts_;
        }

        List<label>& boundaryPatchSizes()
        {
            return boundaryPatchSizes_;
        }

        pointField& points()
        {
            return points_;
        }

        faceList& faces()
        {
            return faces_;
        }

        labelList& owner()
        {
            return owner_;
        }

        labelList& neighbour()
        {
            return neighbour_;
        }

        // fields
        scalarField cellLevelField()
        {
            return cellLevelField_;
        }

        scalarField distanceField()
        {
            return distanceField_;
        }

        scalarField patchIDField()
        {
            return patchIDField_;
        }

        scalarField procIDField()
        {
            return procIDField_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const VDBGridsConverter&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
