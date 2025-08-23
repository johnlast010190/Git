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
    (c) 2011-2022 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/polyTopoChangeMap/mapAddedPolyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mapAddedPolyMesh::mapAddedPolyMesh
(
    const label nOldPoints,
    const label nOldFaces,
    const label nOldCells,
    const label nAddedPoints,
    const label nAddedFaces,
    const label nAddedCells,
    const labelList& oldPointMap,
    const labelList& oldFaceMap,
    const labelList& oldCellMap,

    const labelList& addedPointMap,
    const labelList& addedFaceMap,
    const labelList& addedCellMap,

    const labelList& oldPatchMap,
    const labelList& addedPatchMap,
    const labelList& oldPatchSizes,
    const labelList& oldPatchStarts
)
:
    nOldPoints_(nOldPoints),
    nOldFaces_(nOldFaces),
    nOldCells_(nOldCells),
    nAddedPoints_(nAddedPoints),
    nAddedFaces_(nAddedFaces),
    nAddedCells_(nAddedCells),

    oldPointMap_(oldPointMap),
    oldFaceMap_(oldFaceMap),
    oldCellMap_(oldCellMap),

    addedPointMap_(addedPointMap),
    addedFaceMap_(addedFaceMap),
    addedCellMap_(addedCellMap),

    oldPatchMap_(oldPatchMap),
    addedPatchMap_(addedPatchMap),

    oldPatchSizes_(oldPatchSizes),
    oldPatchStarts_(oldPatchStarts)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
