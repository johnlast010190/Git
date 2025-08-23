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

#include "faceZoneDataSetProvider.H"
#include "Utils/vtkDataHandlingTools.H"

#include "vtk/adaptor/foamVtkVtuAdaptor.H"
#include "meshes/polyMesh/zones/faceZone/faceZone.H"
#include "fields/GeometricFields/GeometricField/GeometricField.H"
#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"

#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkWrappers/engysLabelCoreDataArray.H"
#include "Utils/patchTypes.H"

namespace Foam::functionObjects::runTimeVis
{

template<class fieldType>
void FaceZoneDataSetProvider::addFaceZoneDataToConvertedMesh
    (
        const word &foamFieldName,
        const faceZone &myFaceZone,
        vtkDataSet* convertedMesh
    )
{
    // Either a volScalarField or a volVectorField
    const auto &myField = meshAndFields_.getFieldsRegistry().lookupObject<VolField<fieldType>>(
        foamFieldName
    );

    Field < fieldType > subsettedField(myFaceZone.size());
    const fvMesh &mesh = meshAndFields_.getMesh();
    const List<label> &faceNeighbour = mesh.faceNeighbour();
    const List<label> &faceOwner = mesh.faceOwner();
    label nNeighbors = faceNeighbour.size();
    label nInternalFaces = mesh.nInternalFaces();

    for (label i = 0; i < myFaceZone.size(); ++i)
    {
        label faceId = myFaceZone[i];
        if (faceId >= nInternalFaces)
        {
            // this is a boundary cell, take the owner
            label patchIndex = mesh.boundaryMesh().whichPatch(faceId);
            const fvPatch& patch = mesh.boundary()[patchIndex];

            if (PatchTypes::isProcessPatch(patch))
            {
                label patchFaceId = faceId - patch.start();
                label cellId = faceOwner[faceId];
                subsettedField[i] = (myField[cellId] + myField.boundaryField()[patchIndex][patchFaceId]) * 0.5;
            }
            else if (PatchTypes::isEmptyPatch(patch) || PatchTypes::isNCCPatch(patch))
            {
                label cellId = faceOwner[faceId];
                subsettedField[i] = myField[cellId];
            }
            else
            {
                label patchFaceId = faceId - patch.start();
                subsettedField[i] = myField.boundaryField()[patchIndex][patchFaceId];
            }
        }
        else
        {
            label cellId = faceOwner[faceId];
            if (faceId < nNeighbors)
            {
                label neighborId = faceNeighbour[faceId];
                subsettedField[i] = (myField[cellId] + myField[neighborId]) * 0.5;
            }
            else
            {
                subsettedField[i] = myField[cellId];
            }
        }
    }

    vtkSmartPointer<vtkFloatArray> floatArray =
        Foam::vtk::Tools::convertFieldToVTK<fieldType>
            (
                foamFieldName,
                subsettedField
            );

    if (floatArray->GetNumberOfTuples() == convertedMesh->GetNumberOfCells())
    {
        convertedMesh->GetCellData()->AddArray(floatArray);
    }
    else if (floatArray->GetNumberOfTuples() == convertedMesh->GetNumberOfPoints())
    {
        // Only has vertices
        convertedMesh->GetPointData()->AddArray(floatArray);
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void FaceZoneDataSetProvider::update(scalar vtkNotUsed(currentTime))
{
    const faceZoneMesh &myFaceZoneMesh = meshAndFields_.getMesh().faceZones();
    string name(name_);
    const label zoneIndex = myFaceZoneMesh.findIndex(name);

    if (zoneIndex < 0)
    {
        FatalErrorInFunction << "Error reading face zone " << name << ": "
                             << "not found" << abort(FatalError);
    }

    const faceZone &myFaceZone = myFaceZoneMesh[zoneIndex];
    const primitiveFacePatch& facePatch = myFaceZone();

    // Set geometric data
    vtkSmartPointer<vtkDataSet> convertedMesh = Foam::vtk::Tools::Patch::mesh(facePatch);

    for (const foamField &field: getRequiredFields().getFoamFields())
    {
        if (meshAndFields_.getFieldsRegistry().foundObject<volScalarField>(field))
        {
            addFaceZoneDataToConvertedMesh<scalar>(field, myFaceZone, convertedMesh);
        }
        else if (meshAndFields_.getFieldsRegistry().foundObject<volVectorField>(field))
        {
            addFaceZoneDataToConvertedMesh<vector>(field, myFaceZone, convertedMesh);
        }
        else
        {
            addToNonExistingFields(field);
        }
    }

    addInterpolatedCellAndPointData(convertedMesh);

    vtkNew<engysLabelCoreDataArray> meshPointsArray;
    meshPointsArray->SetCorePointer(&facePatch.meshPoints());
    meshPointsArray->SetName("engysGhostMeshPointsArray");
    convertedMesh->GetPointData()->AddArray(meshPointsArray);

    output_ = setGhostCells(convertedMesh);
}

} // End namespace


// ************************************************************************* //
