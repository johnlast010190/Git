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
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::cellCellStencils::cellVolumeWeight

Description
    Volume-weighted interpolation stencil

SourceFiles
    cellVolumeWeightCellCellStencil.C
    cellVolumeWeightCellCellStencilTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef cellCellStencils_cellVolumeWeight_H
#define cellCellStencils_cellVolumeWeight_H

#include "cellCellStencil/cellCellStencil/cellCellStencil.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cellCellStencils
{

/*---------------------------------------------------------------------------*\
                     Class cellVolumeWeight Declaration
\*---------------------------------------------------------------------------*/

class cellVolumeWeight
:
    public cellCellStencil
{
protected:

    // Static data members

        //- Default overlap tolerance. Fraction of volume
        static scalar defaultOverlapTolerance_;


    // Protected data

        //- Dictionary of motion control parameters
        const dictionary dict_;

        //- Tolerance for volume overlap. Fraction of volume
        scalar overlapTolerance_;

        //- Per cell the cell type
        labelList cellTypes_;

        //- Indices of interpolated cells
        labelList interpolationCells_;

        //- Fetch interpolated cells
        distributionMap cellInterpolationMap_;

        //- Interpolation stencil
        labelListList cellStencil_;

        //- Interpolation weights
        List<scalarList> cellInterpolationWeights_;

        //- Amount of interpolation
        volScalarField cellInterpolationWeight_;


   // Protected Member Functions

        void walkFront
        (
            const scalar layerRelax,
            labelList& allCellTypes,
            scalarField& allWeight
        ) const;

        //- Find cells next to cells of type PATCH
        void findHoles
        (
            const globalIndex& globalCells,
            const fvMesh& mesh,
            const labelList& zoneID,
            const labelListList& stencil,
            labelList& cellTypes
        ) const;

        //- according to additionalDocumentation/MEJ_oversetMesh.txt
        void markPatchCells
        (
            const fvMesh& mesh,
            const labelList& cellMap,
            labelList& patchCellTypes
        ) const;

        void combineCellTypes
        (
            const label subZoneID,
            const fvMesh& subMesh,
            const labelList& subCellMap,

            const label donorZoneID,
            const labelListList& toOtherCells,
            const List<scalarList>& weights,
            const labelList& otherCells,
            const labelList& interpolatedOtherPatchTypes,

            labelListList& allStencil,
            scalarListList& allWeights,
            labelList& allCellTypes,
            labelList& allDonorID
        ) const;

        //- interpolate (= combine) patch types
        void interpolatePatchTypes
        (
            const labelListList& addressing,
            const labelList& patchTypes,
            labelList& result
        ) const;

        //- interpolate (= combine) patch types
        void interpolatePatchTypes
        (
            const autoPtr<distributionMap>& mapPtr,
            const labelListList& addressing,
            const labelList& patchTypes,
            labelList& result
        ) const;


public:

    //- Runtime type information
    TypeName("cellVolumeWeight");


    // Constructors

        //- Construct from fvMesh
        cellVolumeWeight(const fvMesh&, const dictionary&, const bool doUpdate);

        //- Disallow default bitwise copy construction
        cellVolumeWeight(const cellVolumeWeight&) = delete;


    //- Destructor
    virtual ~cellVolumeWeight();


    // Member Functions

        //- Access to volume overlap tolerance
        scalar overlapTolerance() const
        {
            return overlapTolerance_;
        }

        //- Update stencils. Return false if nothing changed.
        virtual bool update();

        //- Return the cell type list
        virtual const labelUList& cellTypes() const
        {
            return cellTypes_;
        }

        //- Indices of interpolated cells
        virtual const labelUList& interpolationCells() const
        {
            return interpolationCells_;
        }

        //- Return a communication schedule
        virtual const distributionMap& cellInterpolationMap() const
        {
            return cellInterpolationMap_;
        }

        //- Per interpolated cell the neighbour cells (in terms of slots as
        //  constructed by above cellInterpolationMap) to interpolate
        virtual const labelListList& cellStencil() const
        {
            return cellStencil_;
        }

        //- Weights for cellStencil
        virtual const List<scalarList>& cellInterpolationWeights() const
        {
            return cellInterpolationWeights_;
        }

        //- Per interpolated cell the interpolation factor. (0 = use
        //  calculated, 1 = use interpolated)
        virtual const scalarList& cellInterpolationWeight() const
        {
            return cellInterpolationWeight_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cellVolumeWeight&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cellCellStencils
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
