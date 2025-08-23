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
    (c) 2021 Engys Ltd.

Class
    Foam::autoLayerCellsMerge

Description
    Merge layers stack cells when in error

SourceFiles
    autoLayerCellsMerge.C

\*---------------------------------------------------------------------------*/

#ifndef autoLayerCellsMerge_H
#define autoLayerCellsMerge_H

#include "meshRefinement/meshRefinement.H"
#include "snappyHexMeshDriver/layerParameters/layerParameters.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class autoLayerCellsMerge Declaration
\*---------------------------------------------------------------------------*/

class autoLayerCellsMerge
{
    // Private data

        //- Mesh+level information
        meshRefinement& meshRefiner_;

        //- Reference to decomposition method
        decompositionMethod& decomposer_;

        //- Reference to mesh distribution engine
        meshDistribute& distributor_;

        //- Container for layer information
        const layerParameters& layerParams_;

    // Private Member Functions

public:

    //- Runtime type information
    ClassName("autoLayerCellsMerge");

    //mark layer cells in error and re-distribute
    // so merge stacks on same processor
    bool markAndRedistribute
    (
        boolList& markedCells
    );

    //Perform merging of cells and update the mesh
    label topoChange
    (
        const boolList& markedCells
    );

    //Clean up hanging points after cell merging
    void removePoints();

    //- top level merge
    void merge
    (
        const label nMergeIter
    );

    // Constructors
    autoLayerCellsMerge
    (
        meshRefinement& meshRefiner,
        decompositionMethod& decomposer,
        meshDistribute& distributor,
        const layerParameters& layerParams
    );

    //- Destructor
    ~autoLayerCellsMerge();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
