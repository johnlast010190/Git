/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.3.1
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
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "patchTypes.H"

#include "meshes/polyMesh/polyPatches/constraint/empty/emptyPolyPatch.H"
#include "meshes/polyMesh/polyPatches/derived/nonConformalOrig/nonConformalOrigPolyPatch.H"
#include "nonConformal/polyPatches/nonConformalError/nonConformalErrorPolyPatch.H"
#include "nonConformal/polyPatches/nonConformalProcessorCyclic/nonConformalProcessorCyclicPolyPatch.H"
#include "fvMesh/fvPatches/constraint/internal/internalFvPatch.H"
#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"

namespace Foam::functionObjects::runTimeVis
{

bool PatchTypes::isEmptyPatch(const fvPatch& patch)
{
    return isType<emptyPolyPatch>(patch.patch());
}

bool PatchTypes::isCoupledPatch(const fvPatch& patch)
{
    return isType<nonConformalOrigPolyPatch>(patch.patch()) || patch.coupled();
}

bool PatchTypes::isNCCPatch(const fvPatch& patch)
{
    const polyPatch& polyPath = patch.patch();
    return isType<nonConformalOrigPolyPatch>(polyPath) ||
           isType<nonConformalCyclicPolyPatch>(polyPath) ||
           isType<nonConformalErrorPolyPatch>(polyPath) ||
           isType<nonConformalProcessorCyclicPolyPatch>(polyPath);
}

bool PatchTypes::isProcessPatch(const fvPatch& patch)
{
    return patch.patch().type() == processorPolyPatch::typeName_();
}

bool PatchTypes::isInternalPatch(const fvPatch &patch)
{
    return isA<internalFvPatch>(patch) || isA<indirectPrimitivePatch>(patch.patch());
}

} // End namespace
