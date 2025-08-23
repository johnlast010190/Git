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

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

\*---------------------------------------------------------------------------*/

#include "utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressing.H"
#include "include/demandDrivenData.H"
#include "meshes/primitiveMesh/primitiveMesh.H"

# ifdef USE_OMP
#include <omp.h>
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyMeshGenAddressing::updateGeometry
(
    const boolList& changedFace
)
{
    const pointFieldPMG& p = mesh_.points();
    const faceListPMG& faces = mesh_.faces();

    //- update face centres and face areas
    if (faceCentresPtr_ && faceAreasPtr_)
    {
        vectorField& fCtrs = *faceCentresPtr_;
        vectorField& fAreas = *faceAreasPtr_;

        # ifdef USE_OMP
        # pragma omp parallel for if (faces.size() > 100) \
        schedule(dynamic, 10)
        # endif
        forAll(faces, facei)
        {
            if (changedFace[facei])
            {
                const face& f = faces[facei];
                primitiveMesh::faceAreaAndCentre(f,p,fCtrs[facei],fAreas[facei]);
            }
        }
    }

    //- update cell centres and cell volumes
    if (cellCentresPtr_ && cellVolumesPtr_ && faceCentresPtr_ && faceAreasPtr_)
    {
        const vectorField& fCtrs = *faceCentresPtr_;
        const vectorField& fAreas = *faceAreasPtr_;
        vectorField& cellCtrs = *cellCentresPtr_;
        scalarField& cellVols = *cellVolumesPtr_;

        const labelList& own = mesh_.owner();
        const cellListPMG& cells = mesh_.cells();

        # ifdef USE_OMP
        # pragma omp parallel for if (cells.size() > 100) \
        schedule(dynamic, 10)
        # endif
        forAll(cells, celli)
        {
            const cell& c = cells[celli];

            bool update(false);
            forAll(c, fI)
                if (changedFace[c[fI]])
                {
                    update = true;
                    break;
                }

            if (update)
            {
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
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
