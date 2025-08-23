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

Description
    Calulate the face centres and areas.

    Calculate the centre by breaking the face into triangles using the face
    centre and area-weighted averaging their centres.  This method copes with
    small face-concavity.

\*---------------------------------------------------------------------------*/

#include "utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressing.H"
#include "meshes/primitiveMesh/primitiveMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyMeshGenAddressing::calcFaceCentresAndAreas() const
{
    if (faceCentresPtr_ || faceAreasPtr_)
    {
        FatalErrorInFunction
            << "Face centres or face areas already calculated"
            << abort(FatalError);
    }

    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();

    faceCentresPtr_ = new vectorField(faces.size());
    vectorField& fCtrs = *faceCentresPtr_;

    faceAreasPtr_ = new vectorField(faces.size());
    vectorField& fAreas = *faceAreasPtr_;

    makeFaceCentresAndAreas(points, fCtrs, fAreas);
}

void polyMeshGenAddressing::makeFaceCentresAndAreas
(
    const pointFieldPMG& p,
    vectorField& fCtrs,
    vectorField& fAreas
) const
{
    const faceListPMG& fs = mesh_.faces();
    const label nFaces = fs.size();

    # ifdef USE_OMP
    # pragma omp parallel for if (nFaces > 1000)
    # endif
    for (label facei=0;facei<nFaces;++facei)
    {
        face f = fs[facei];
        primitiveMesh::faceAreaAndCentre(f,p,fCtrs[facei],fAreas[facei]);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const vectorField& polyMeshGenAddressing::faceCentres() const
{
    if (!faceCentresPtr_)
    {
        # ifdef USE_OMP
        if (omp_in_parallel())
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcFaceCentresAndAreas();
    }

    return *faceCentresPtr_;
}

const vectorField& polyMeshGenAddressing::faceAreas() const
{
    if (!faceAreasPtr_)
    {
        # ifdef USE_OMP
        if (omp_in_parallel())
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcFaceCentresAndAreas();
    }

    return *faceAreasPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
