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

#include "volMeshDataSetProvider.H"

#include "Utils/Utils.H"
#include "Utils/basicProfiler.H"
#include "vtkWrappers/rtppVtkVtuAdaptor.H"

#include "unordered_map"

#include "fields/fvPatchFields/fvPatchField/fvPatchField.H"

#define PROFILE_POINT(LABEL)
//#define PROFILE_POINT(label) basicProfiler::addMeasurePoint(label)

namespace Foam::functionObjects::runTimeVis
{

VolMeshDataSetProvider::VolMeshDataSetProvider(
    const string &name,
    const string &region,
    const MeshAndFields &meshAndFields
)
        : MeshDataSetProvider(name, region, meshAndFields)
{
}

VolMeshDataSetProvider::~VolMeshDataSetProvider() = default;

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

class CellToPointFieldBoundaryOverwriter
{
protected:
    const std::vector<const polyPatch*>& boundaries;

    CellToPointFieldBoundaryOverwriter(
        const std::vector<const polyPatch *> &boundaries,
        const std::unordered_map<label, label> &mappedWeights
    ) : boundaries(boundaries), totalWeight(mappedWeights.size()), pointIndexes(mappedWeights.size())
    {
        label i = 0;
        for (const auto& mappedWeight : mappedWeights)
        {
            pointIndexes[i] = mappedWeight.first;
            totalWeight[i] = mappedWeight.second;
            i++;
        }
    };

    template<int nComponents>
    void initializeField(vtkFloatArray* pointArray)
    {
        float* pointArrayPtr = pointArray->GetPointer(0);
        for (label pointId : pointIndexes)
        {
            const label index = pointId * nComponents;
            for (int component = 0; component < nComponents; component++)
            {
                pointArrayPtr[index + component] = 0;
            }
        }
    }

    template<int nComponents>
    void finalizeField(vtkFloatArray* pointArray)
    {
        float* pointArrayPtr = pointArray->GetPointer(0);
        const auto nIndexedPoints = static_cast<label>(pointIndexes.size());
        for (label i = 0; i < nIndexedPoints; i++)
        {
            const label pointId = pointIndexes[i];
            const auto pointWeight = static_cast<float>(totalWeight[i]);
            const label index = pointId * nComponents;
            for (int component = 0; component < nComponents; component++)
            {
                pointArrayPtr[index + component] /= pointWeight;
            }
        }
    }

private:
    std::vector<label> totalWeight;
    std::vector<label> pointIndexes;
};

class CellToPointFieldExternalBoundaryOverwriter : public CellToPointFieldBoundaryOverwriter
{
public:
    explicit CellToPointFieldExternalBoundaryOverwriter(const std::vector<const polyPatch *> &boundaries)
        : CellToPointFieldBoundaryOverwriter(
        boundaries,
        mapPointWeights(boundaries))
    {};

    void overwriteBoundaryValues(vtkFloatArray* pointArray, const std::vector<vtkDataArray*>& vtkPointFields)
    {
        switch(pointArray->GetNumberOfComponents())
        {
            case 1:
                initializeField<1>(pointArray);
                addBoundaryValues<1>(pointArray, vtkPointFields);
                finalizeField<1>(pointArray);
                break;
            case 3:
                initializeField<3>(pointArray);
                addBoundaryValues<3>(pointArray, vtkPointFields);
                finalizeField<3>(pointArray);
                break;
        }
    }

private:
    static std::unordered_map<label, label> mapPointWeights(const std::vector<const polyPatch*>& boundaries)
    {
        label totalBoundarySize = 0;
        for (const polyPatch* boundary : boundaries)
        {
            totalBoundarySize += boundary->size();
        }

        std::unordered_map<label, label> mappedWeights(totalBoundarySize);
        for (const polyPatch* boundary : boundaries)
        {
            const List<label>& meshPoints = boundary->meshPoints();
            for (label pointId : meshPoints)
            {
                mappedWeights[pointId]++;
            }
        }
        return mappedWeights;
    }

    template<int nComponents>
    void addBoundaryValues(vtkFloatArray* pointArray, const std::vector<vtkDataArray*>& vtkPointFields)
    {
        float* pointArrayPtr = pointArray->GetPointer(0);
        label boundaryIndex = 0;
        for (const polyPatch* boundary : boundaries)
        {
            vtkFloatArray* vtkBoundaryPointField = vtkFloatArray::SafeDownCast(vtkPointFields[boundaryIndex]);
            if (!vtkBoundaryPointField) continue;
            float* vtkBoundaryPointFieldPtr = vtkBoundaryPointField->GetPointer(0);
            boundaryIndex++;
            const polyPatch& patch = *boundary;
            const UList<label>& meshPoints = patch.meshPoints();
            const label nBoundaryPoints = meshPoints.size();
            for (label pointIndex = 0; pointIndex < nBoundaryPoints; pointIndex++)
            {
                const label pointId = meshPoints[pointIndex];
                const label internalMeshIndex = pointId * nComponents;
                const label boundaryMeshIndex = pointIndex * nComponents;
                for (int component = 0; component < nComponents; component++)
                {
                    pointArrayPtr[internalMeshIndex + component] += vtkBoundaryPointFieldPtr[boundaryMeshIndex + component];
                }
            }
        }
    }
};

void VolMeshDataSetProvider::update(scalar vtkNotUsed(currentTime))
{
    // Calling vtuAdaptor::internal constructs an unstructured grid from the foam
    // mesh.
    PROFILE_POINT("updating volmesh");
    const fvMesh& coreMesh = meshAndFields_.getMesh();
    rtppVtuAdaptor rtppVtuAdaptor(adaptor_);
    vtkSmartPointer<vtkUnstructuredGrid> convertedMesh = rtppVtuAdaptor.internal(coreMesh);

    if (!convertedMesh)
    {
        FatalError << "Failure converting the FOAM volume mesh to vtk: "
            << name_
            << exit(FatalError);
    }
    PROFILE_POINT("converted volmesh");

    for (const foamField& field : getRequiredFields().getFoamFields())
    {
        addField(field, convertedMesh);
    }
    PROFILE_POINT("added fields");

    addInterpolatedCellAndPointData(convertedMesh);
    PROFILE_POINT("interpolated points");

    if (getRequiredFields().hasPointField())
    {
        // This overwrites the values for points in the external boundaries by interpolating only from boundary cells
        overridePointValuesWithBoundaryPointValues(convertedMesh);
        PROFILE_POINT("overrode points");
    }

    if (!convertedMesh)
    {
        FatalError << "Interpolation of point values failed for volume mesh "
                   << name_
                   << exit(FatalError);
    }

    output_ = setGhostCells(convertedMesh, true);
    PROFILE_POINT("set ghost cells");
}

bool VolMeshDataSetProvider::addField
(
        const foamField& fieldName,
        vtkDataSet* convertedMesh
)
{
    bool isPointArray;
    vtkSmartPointer<vtkDataArray> vtkData = getVolmeshFieldDataArray(fieldName, isPointArray);

    if (vtkData)
    {
        if (!isPointArray && vtkData->GetNumberOfTuples() == convertedMesh->GetNumberOfCells())
        {
            convertedMesh->GetCellData()->AddArray(vtkData);
            return true;
        }
        else if (isPointArray && vtkData->GetNumberOfTuples() == convertedMesh->GetNumberOfPoints())
        {
            convertedMesh->GetPointData()->AddArray(vtkData);
            nativePointFields_.insert(vtkData->GetName());
            return true;
        }
    }
    return false;
}

void VolMeshDataSetProvider::overridePointValuesWithBoundaryPointValues(vtkDataSet *internalDataSet)
{
    const fvMesh& coreMesh = meshAndFields_.getMesh();

    std::vector<const polyPatch*> externalBoundaries;
    std::vector<ItemDataSetProvider*> externalBoundaryProviders;
    for (ItemDataSetProvider* source : getUsedSources())
    {
        label patchId = meshAndFields_.getMesh().boundary().findPatchID(source->getName());

        if (patchId < 0) continue;
        const polyPatch &boundary = coreMesh.boundaryMesh()[patchId];
        if (boundary.meshPoints().size() <= 0) continue;

        if (source->isProcessorBoundary())
        {
            continue;
        }
        else if (!source->isNCCPatch() &&
                 !source->isCoupledPatch() &&
                 !source->isEmptyPatch())
        {
            externalBoundaryProviders.push_back(source);
            externalBoundaries.push_back(&boundary);
        }
    }

    CellToPointFieldExternalBoundaryOverwriter externalBoundaryOverwriter(externalBoundaries);
    const auto nExternalBoundaryProviders = static_cast<label>(externalBoundaryProviders.size());

    std::vector<vtkDataArray*> externalBoundaryVtkFields(nExternalBoundaryProviders);

    for (int arrayId = 0; arrayId < internalDataSet->GetCellData()->GetNumberOfArrays(); arrayId++)
    {
        vtkDataArray* cellArray = internalDataSet->GetCellData()->GetArray(arrayId);
        const char* arrayName = cellArray->GetName();
        if (nativePointFields_.find(arrayName) != nativePointFields_.end()) continue;
        vtkFloatArray* pointArray = vtkFloatArray::SafeDownCast(internalDataSet->GetPointData()->GetArray(arrayName));
        if (!pointArray) continue;

        for (label externalBoundaryIndex = 0; externalBoundaryIndex < nExternalBoundaryProviders; externalBoundaryIndex++)
        {
            externalBoundaryVtkFields[externalBoundaryIndex] = externalBoundaryProviders[externalBoundaryIndex]->getDataSetOutput()->GetPointData()->GetArray(arrayName);
        }

        externalBoundaryOverwriter.overwriteBoundaryValues(pointArray, externalBoundaryVtkFields);
    }
}

} // End namespace


// ************************************************************************* //
