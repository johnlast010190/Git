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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fvMesh/extendedStencil/faceToCell/extendedCentredFaceToCellStencil.H"
#include "meshes/polyMesh/polyDistributionMap/distributionMap.H"
#include "fvMesh/extendedStencil/faceToCell/globalIndexStencils/faceToCellStencil.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedCentredFaceToCellStencil::extendedCentredFaceToCellStencil
(
    const faceToCellStencil& stencil
)
:
    extendedFaceToCellStencil(stencil.mesh()),
    stencil_(stencil)
{
    // Calculate distribute map (also renumbers elements in stencil)
    List<Map<label>> compactMap(Pstream::nProcs());
    mapPtr_.reset
    (
        new distributionMap
        (
            stencil.globalNumbering(),
            stencil_,
            compactMap
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::extendedCentredFaceToCellStencil::compact()
{
    boolList isInStencil(map().constructSize(), false);

    forAll(stencil_, facei)
    {
        const labelList& stencilCells = stencil_[facei];

        forAll(stencilCells, i)
        {
            isInStencil[stencilCells[i]] = true;
        }
    }

    mapPtr_().compact(isInStencil, Pstream::msgType());
}


// ************************************************************************* //
