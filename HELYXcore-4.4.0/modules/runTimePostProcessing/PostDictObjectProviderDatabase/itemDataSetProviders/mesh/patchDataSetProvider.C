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
    (c) 2015-2019 OpenCFD Ltd.
    (c) 2022-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "patchDataSetProvider.H"

#include "postDict/postDictKeys.H"
#include "Utils/patchTypes.H"

#include "vtkBitArray.h"
#include "vtkMultiProcessController.h"
#include "vtkWrappers/engysLabelCoreDataArray.H"
#include "meshes/polyMesh/polyPatches/indirectPolyPatch/indirectPolyPatch.H"

namespace Foam::functionObjects::runTimeVis
{

PatchDataSetProvider::PatchDataSetProvider(const string &name, const string &region, const MeshAndFields& meshAndFields)
    : MeshDataSetProvider(name, region, meshAndFields)
{
}

PatchDataSetProvider::~PatchDataSetProvider() = default;

// * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * * //

vtkSmartPointer<vtkFloatArray> PatchDataSetProvider::getConvertedFieldForPatch
    (
        const foamField &foamFieldName,
        const label &patchIndex,
        bool& isPointArray
    )
{
    vtkSmartPointer<vtkFloatArray> field = meshAndFields_.convertBoundaryFoamFieldToVtkArray(foamFieldName,
        patchIndex,
        isPointArray
    );
    if (!field)
    {
        addToNonExistingFields(foamFieldName);
    }
    return field;
}

vtkSmartPointer<vtkFloatArray> PatchDataSetProvider::getConvertedFieldForEmptyPatch
    (
        const foamField &foamFieldName,
        const fvPatch& patch,
        bool& isPointArray
    )
{
    vtkSmartPointer<vtkFloatArray> field = meshAndFields_.convertEmptyBoundaryFoamFieldToVtkArray(foamFieldName,
                                                                                             patch,
                                                                                             isPointArray
    );
    if (!field)
    {
        addToNonExistingFields(foamFieldName);
    }
    return field;
}

void PatchDataSetProvider::addArrayToConvertedMesh
    (
        const foamField &fieldName,
        const fvPatch& patch,
        const label &patchIndex,
        vtkDataSet* convertedMesh
    )
{
    vtkSmartPointer<vtkFloatArray> floatArray;
    bool isPointArray;
    if (isEmptyPatch() || isNCCPatch())
    {
        floatArray = getConvertedFieldForEmptyPatch(fieldName, patch, isPointArray);
    }
    else
    {
        floatArray = getConvertedFieldForPatch(fieldName, patchIndex, isPointArray);
    }
    if (floatArray)
    {
        if (isPointArray && floatArray->GetNumberOfTuples() == convertedMesh->GetNumberOfPoints())
        {
            convertedMesh->GetPointData()->AddArray(floatArray);
        }
        else if (floatArray->GetNumberOfTuples() == convertedMesh->GetNumberOfCells())
        {
            convertedMesh->GetCellData()->AddArray(floatArray);
        }
    }
}

void PatchDataSetProvider::update(scalar vtkNotUsed(currentTime))
{
    label patchIndex = meshAndFields_.getMesh().boundary().findPatchID(name_);
    if (patchIndex < 0)
    {
        // Patch not found in this processor
        vtkNew<vtkPolyData> emptyData;
        output_ = setGhostCells(emptyData);
        return;
    }
    const fvPatch &foamPatch = meshAndFields_.getMesh().boundary()[patchIndex];
    const polyPatch& foamPolyPatch = foamPatch.patch();
    isEmptyPatch_ = PatchTypes::isEmptyPatch(foamPatch);
    isCoupledPatch_ = PatchTypes::isCoupledPatch(foamPatch);
    isNCCPatch_ = PatchTypes::isNCCPatch(foamPatch);

    vtkSmartPointer<vtkPolyData> convertedPatch = Foam::vtk::Tools::Patch::mesh(foamPolyPatch);

    if (!convertedPatch)
    {
        Warning << "Error converting patch " << name_ << " to vtk" << endl;
    }

    for (const foamField &field: getRequiredFields().getFoamFields())
    {
        // Add the colour field to the patch
        addArrayToConvertedMesh(foamField(field), foamPatch, patchIndex, convertedPatch);
    }

    addInterpolatedCellAndPointData(convertedPatch);

    vtkNew<engysLabelCoreDataArray> meshPointsArray;
    meshPointsArray->SetCorePointer(&foamPolyPatch.meshPoints());
    meshPointsArray->SetName("engysGhostMeshPointsArray");
    convertedPatch->GetPointData()->AddArray(meshPointsArray);

    output_ = setGhostCells(convertedPatch);

    if (isA<indirectPolyPatch>(foamPolyPatch))
    {
        handleConditionalPatch(foamPatch);
    }
}

void PatchDataSetProvider::handleConditionalPatch(const fvPatch &conditionalPatch)
{
    bool active = conditionalPatch.isActive();

    vtkSmartPointer<vtkDataArray> bitArray = output_->GetFieldData()->GetArray(activePatchKeys::ACTIVE_FIELD_KEY);
    if (bitArray == nullptr)
    {
        bitArray = vtkSmartPointer<vtkBitArray>::New();
        bitArray->SetNumberOfComponents(1);
        bitArray->Allocate(1);
        bitArray->SetName(activePatchKeys::ACTIVE_FIELD_KEY);
        output_->GetFieldData()->AddArray(bitArray);
    }
    bitArray->InsertTuple1(0, active);
}

} // End namespace


// ************************************************************************* //
