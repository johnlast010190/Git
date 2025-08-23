/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.0.1
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
    (c) 2020-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "meshAndFields.H"

#include "foamField.H"
#include "vtk/adaptor/foamVtkVtuAdaptor.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchField.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchField.H"
#include "meshes/pointMesh/pointMesh.H"
#include "fields/GeometricFields/pointFields/pointFieldsFwd.H"

#include "vtkWrappers/engysScalarCoreDataArray.H"
#include "vtkWrappers/engysVectorCoreDataArray.H"

namespace Foam::functionObjects::runTimeVis
{

MeshAndFields::VectorOrScalar MeshAndFields::getFieldType(const std::string &foamFieldName) const
{
    if (fieldsRegistry.foundObject<volScalarField>(foamFieldName))
    {
        return SCALAR;
    }
    else if (fieldsRegistry.foundObject<volVectorField>(foamFieldName))
    {
        return VECTOR;
    }
    else if (fieldsRegistry.foundObject<pointScalarField>(foamFieldName))
    {
        return POINT_SCALAR;
    }
    else if (fieldsRegistry.foundObject<pointVectorField>(foamFieldName))
    {
        return POINT_VECTOR;
    }
    else
    {
        return NOT_FOUND;
    }
}

bool MeshAndFields::isMeshFieldValid(const foamField &fieldName) const
{
    return getFieldType(fieldName.getFoamName()) != NOT_FOUND;
}

template<class type>
vtkSmartPointer<vtkDataArray> MeshAndFields::convertInternal(const foamField &fieldName) const
{
    const auto &field = fieldsRegistry.lookupObject<type>(fieldName.getFoamName());
    return vtk::Tools::convertFieldToVTK(fieldName, field);
}

vtkSmartPointer<vtkDataArray> MeshAndFields::convertInternalScalar(const foamField &fieldName) const
{
    const auto &field = fieldsRegistry.lookupObject<volScalarField>(fieldName.getFoamName());
    vtkNew<engysScalarCoreDataArray> wrappedArray;
    wrappedArray->SetCorePointer(&field);
    wrappedArray->SetName(fieldName.getFoamName().c_str());
    return wrappedArray;
}

vtkSmartPointer<vtkDataArray> MeshAndFields::convertInternalVector(const foamField &fieldName) const
{
    const auto &field = fieldsRegistry.lookupObject<volVectorField>(fieldName.getFoamName());
    vtkNew<engysVectorCoreDataArray> wrappedArray;
    wrappedArray->SetCorePointer(&field);
    wrappedArray->SetName(fieldName.getFoamName().c_str());
    return wrappedArray;
}

vtkSmartPointer<vtkDataArray> MeshAndFields::convertInternalFoamFieldToVtkArray(
    const foamField &fieldName,
    bool& isPointArray) const
{
    const std::string &foamFieldName = fieldName.getFoamName();
    switch (getFieldType(foamFieldName))
    {
        case SCALAR:
            isPointArray = false;
            return convertInternalScalar(fieldName);
        case POINT_SCALAR:
            isPointArray = true;
            return convertInternal<pointScalarField>(fieldName);
        case VECTOR:
            isPointArray = false;
            return convertInternalVector(fieldName);
        case POINT_VECTOR:
            isPointArray = true;
            return convertInternal<pointVectorField>(fieldName);
        default:
            isPointArray = false;
            return nullptr;
    }
}

template<class type>
vtkSmartPointer<vtkFloatArray> MeshAndFields::convertVolBoundary(const foamField &fieldName, label boundaryIndex) const
{
    const auto &field = fieldsRegistry.lookupObject<type>(fieldName.getFoamName()).boundaryField()[boundaryIndex];
    return vtk::Tools::convertFieldToVTK(fieldName, field);
}

template<class type>
vtkSmartPointer<vtkFloatArray> MeshAndFields::convertPointBoundary(const foamField &fieldName, label boundaryIndex) const
{
    const auto &field = fieldsRegistry.lookupObject<type>(fieldName.getFoamName()).boundaryField()[boundaryIndex];
    return vtk::Tools::convertFieldToVTK(fieldName, field.internalField());
}

template<class type, label nComponents>
vtkSmartPointer<vtkFloatArray> MeshAndFields::convertEmptyVolBoundary(const foamField &fieldName, const fvPatch& boundaryPatch) const
{
    const auto &field = fieldsRegistry.lookupObject<type>(fieldName.getFoamName());
    vtkNew<vtkFloatArray> array;
    array->SetName(fieldName.c_str());
    array->SetNumberOfComponents(nComponents);
    label nFaces = boundaryPatch.patch().size();
    array->SetNumberOfTuples(nFaces);

    const List<label>& faceOwner = boundaryPatch.boundaryMesh().mesh().faceOwner();
    const label offset = boundaryPatch.start();

    float* arrayPointer = array->GetPointer(0);

    for (label i = 0; i < nFaces; i++)
    {
        for (label c = 0; c < nComponents; c++)
        {
            *arrayPointer = component(field[faceOwner[i + offset]], c);
            arrayPointer++;
        }
    }

    return array;
}

template<class type>
vtkSmartPointer<vtkFloatArray> MeshAndFields::convertEmptyPointBoundary(const foamField &fieldName, const fvPatch& boundaryPatch) const
{
    WarningInFunction << "Point array for empty patches is not supported" << endl;
    return nullptr;
}

vtkSmartPointer<vtkFloatArray> MeshAndFields::convertBoundaryFoamFieldToVtkArray(
    const foamField &fieldName,
    label boundaryIndex,
    bool& isPointArray
) const
{
    const std::string &foamFieldName = fieldName.getFoamName();
    switch (getFieldType(foamFieldName))
    {
        case SCALAR:
            isPointArray = false;
            return convertVolBoundary<volScalarField>(fieldName, boundaryIndex);
        case POINT_SCALAR:
            isPointArray = true;
            return convertPointBoundary<pointScalarField>(fieldName, boundaryIndex);
        case VECTOR:
            isPointArray = false;
            return convertVolBoundary<volVectorField>(fieldName, boundaryIndex);
        case POINT_VECTOR:
            isPointArray = true;
            return convertPointBoundary<pointVectorField>(fieldName, boundaryIndex);
        default:
            isPointArray = false;
            return nullptr;
    }
}

vtkSmartPointer<vtkFloatArray> MeshAndFields::convertEmptyBoundaryFoamFieldToVtkArray(
    const foamField &fieldName,
    const fvPatch& boundaryPatch,
    bool& isPointArray
) const
{
    const std::string &foamFieldName = fieldName.getFoamName();
    switch (getFieldType(foamFieldName))
    {
        case SCALAR:
            isPointArray = false;
            return convertEmptyVolBoundary<volScalarField, 1>(fieldName, boundaryPatch);
        case POINT_SCALAR:
            isPointArray = true;
            return convertEmptyPointBoundary<pointScalarField>(fieldName, boundaryPatch);
        case VECTOR:
            isPointArray = false;
            return convertEmptyVolBoundary<volVectorField, 3>(fieldName, boundaryPatch);
        case POINT_VECTOR:
            isPointArray = true;
            return convertEmptyPointBoundary<pointVectorField>(fieldName, boundaryPatch);
        default:
            isPointArray = false;
            return nullptr;
    }
}

template<class T>
scalarMinMax MeshAndFields::calculateVectorMinMax(const foamField &fieldName, const T &myField) const
{
    if (fieldName.isComponent())
    {
        // Components
        return {
            min(myField).value()[fieldName.getComponentIndex()],
            max(myField).value()[fieldName.getComponentIndex()]
        };
    }
    else
    {
        // Magnitudes
        return {
            min(mag(myField)).value(),
            max(mag(myField)).value()
        };
    }
}

scalarMinMax MeshAndFields::calculateDomainRangeForVectorField(const foamField &fieldName) const
{
    const auto &myField = fieldsRegistry.lookupObject<volVectorField>(fieldName.getFoamName());
    return calculateVectorMinMax(fieldName, myField);
}

scalarMinMax MeshAndFields::calculateDomainRangeForPointVectorField(const foamField &fieldName) const
{
    const auto &myField = fieldsRegistry.lookupObject<pointVectorField>(fieldName.getFoamName());
    return calculateVectorMinMax(fieldName, myField.internalField());
}

scalarMinMax MeshAndFields::calculateDomainRangeForScalarField(const foamField &fieldName) const
{
    const auto &myField = fieldsRegistry.lookupObject<volScalarField>(fieldName.getFoamName());
    return {
        min(myField).value(),
        max(myField).value()
    };
}

scalarMinMax MeshAndFields::calculateDomainRangeForPointScalarField(const foamField &fieldName) const
{
    const auto &myField = fieldsRegistry.lookupObject<pointScalarField>(fieldName.getFoamName());
    return {
        min(myField.internalField()).value(),
        max(myField.internalField()).value()
    };
}

scalarMinMax MeshAndFields::calculateDomainRangeForField(const foamField &fieldName) const
{
    bool fieldNameIsAVector = fieldName.isComponent() || fieldName.isMagnitude();
    VectorOrScalar actualType = getFieldType(fieldName.getFoamName());

    switch (actualType)
    {
        case SCALAR:
            if (!fieldNameIsAVector)
            {
                return calculateDomainRangeForScalarField(fieldName);
            }
            break;
        case VECTOR:
            if (fieldNameIsAVector)
            {
                return calculateDomainRangeForVectorField(fieldName);
            }
            break;
        case POINT_VECTOR:
            if (fieldNameIsAVector)
            {
                return calculateDomainRangeForPointVectorField(fieldName);
            }
            break;
        case POINT_SCALAR:
            if (!fieldNameIsAVector)
            {
                return calculateDomainRangeForPointScalarField(fieldName);
            }
            break;
        case NOT_FOUND:
            break;
    }
    return {};
}

FoamFields MeshAndFields::listAllFields() const
{
    FoamFields allFields;

    List<word> scalarFields = fieldsRegistry.names<volScalarField>();
    for (const word &field: scalarFields)
    {
        if (shouldFieldBeListed(field))
        {
            allFields.addFoamField(foamField(field));
        }
    }
    List<word> vectorFields = fieldsRegistry.names<volVectorField>();
    for (const word &field: vectorFields)
    {
        if (shouldFieldBeListed(field))
        {
            allFields.addFoamField(foamField(field+"-X"));
            allFields.addFoamField(foamField(field+"-Y"));
            allFields.addFoamField(foamField(field+"-Z"));
            allFields.addFoamField(foamField(field+"-Mag"));
        }
    }
    return allFields;
}

bool MeshAndFields::shouldFieldBeListed(const word &field)
{
    return !field.endsWith("PrevIter") && field.find("_copyBySwak4Foam") == word::npos;
}

}
