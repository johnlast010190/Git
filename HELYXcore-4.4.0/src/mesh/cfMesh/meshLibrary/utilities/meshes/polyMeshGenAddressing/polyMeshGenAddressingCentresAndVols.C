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

Class
    polyMeshGenAddressing

Description
    Efficient cell-centre calculation using face-addressing, face-centres and
    face-areas.

\*---------------------------------------------------------------------------*/

#include "utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressing.H"
#include "meshes/primitiveMesh/primitiveMesh.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyMeshGenAddressing::calcCellCentresAndVols() const
{
    if (cellCentresPtr_ || cellVolumesPtr_)
    {
        FatalErrorInFunction
            << "Cell centres or cell volumes already calculated"
            << abort(FatalError);
    }

    const cellListPMG& cells = mesh_.cells();

    // set the accumulated cell centre to zero vector
    cellCentresPtr_ = new vectorField(cells.size());
    vectorField& cellCtrs = *cellCentresPtr_;

    // Initialise cell volumes to 0
    cellVolumesPtr_ = new scalarField(cells.size());
    scalarField& cellVols = *cellVolumesPtr_;

    // Make centres and volumes
    makeCellCentresAndVols(faceCentres(), faceAreas(), cellCtrs, cellVols);
}

void polyMeshGenAddressing::makeCellCentresAndVols
(
    const vectorField& fCtrs,
    const vectorField& fAreas,
    vectorField& cellCtrs,
    scalarField& cellVols
) const
{
    const labelList& own = mesh_.owner();
    const cellListPMG& cells = mesh_.cells();
    const label nCells = cells.size();

    # ifdef USE_OMP
    # pragma omp parallel for if (nCells > 1000)
    # endif
    for (label celli=0;celli<nCells;++celli)
    {
        const cell& c = cells[celli];

        vector& cellCentre =  cellCtrs[celli];
        scalar& cellVol =  cellVols[celli];

        primitiveMesh::cellCentreAndVol
        (
            celli,
            c,
            own,
            fCtrs,
            fAreas,
            cellCentre,
            cellVol
        );
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const vectorField& polyMeshGenAddressing::cellCentres() const
{
    if (!cellCentresPtr_)
    {
        # ifdef USE_OMP
        if (omp_in_parallel())
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcCellCentresAndVols();
    }

    return *cellCentresPtr_;
}

const scalarField& polyMeshGenAddressing::cellVolumes() const
{
    if (!cellVolumesPtr_)
    {
        # ifdef USE_OMP
        if (omp_in_parallel())
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcCellCentresAndVols();
    }

    return *cellVolumesPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
