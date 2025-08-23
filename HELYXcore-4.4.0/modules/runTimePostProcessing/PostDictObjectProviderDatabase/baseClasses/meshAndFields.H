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
    (c) 2019 OpenCFD Ltd.
    (c) 2023-2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::MeshAndFields

Description
    Struct that holds the fvMesh for the foam geometry and an object registry for the fields

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "primitives/ranges/MinMax/MinMax.H"
#include "baseClasses/foamFields.H"
#include "baseClasses/foamField.H"

#include "vtkSmartPointer.h"

#include <string>

class vtkFloatArray;
class vtkDataArray;

namespace Foam
{
class fvMesh;
class objectRegistry;
class fvPatch;

namespace functionObjects::runTimeVis
{
class foamField;

class MeshAndFields {
public:
    MeshAndFields(const fvMesh &mesh, const objectRegistry &fieldsRegistry)
        : mesh(mesh), fieldsRegistry(fieldsRegistry)
    {};

    MeshAndFields(const MeshAndFields &meshAndFields) = default;

    [[nodiscard]] bool isMeshFieldValid(const foamField &fieldName) const;

    [[nodiscard]] vtkSmartPointer<vtkDataArray> convertInternalFoamFieldToVtkArray(
        const foamField &fieldName,
        bool &isPointArray
    ) const;

    [[nodiscard]] vtkSmartPointer<vtkFloatArray> convertBoundaryFoamFieldToVtkArray(
        const foamField &fieldName,
        label boundaryIndex,
        bool &isPointArray
    ) const;

    [[nodiscard]] vtkSmartPointer<vtkFloatArray> convertEmptyBoundaryFoamFieldToVtkArray(
        const foamField &fieldName,
        const fvPatch& boundaryPatch,
        bool &isPointArray
    ) const;

    [[nodiscard]] scalarMinMax calculateDomainRangeForField(const foamField &fieldName) const;

    [[nodiscard]] FoamFields listAllFields() const;

    [[nodiscard]] const fvMesh& getMesh() const {return mesh;};
    [[nodiscard]] const objectRegistry& getFieldsRegistry() const {return fieldsRegistry;};

private:
    const fvMesh& mesh;
    const objectRegistry& fieldsRegistry;

    enum VectorOrScalar
    {
        VECTOR,
        POINT_VECTOR,
        SCALAR,
        POINT_SCALAR,
        NOT_FOUND
    };

    [[nodiscard]] VectorOrScalar getFieldType(const std::string &foamFieldName) const;

    template<class type>
    [[nodiscard]] vtkSmartPointer<vtkDataArray> convertInternal(const foamField &fieldName) const;
    [[nodiscard]] vtkSmartPointer<vtkDataArray> convertInternalScalar(const foamField &fieldName) const;
    [[nodiscard]] vtkSmartPointer<vtkDataArray> convertInternalVector(const foamField &fieldName) const;

    template<class type>
    [[nodiscard]] vtkSmartPointer<vtkFloatArray> convertVolBoundary(const foamField &fieldName, label boundaryIndex) const;

    template<class type>
    [[nodiscard]] vtkSmartPointer<vtkFloatArray> convertPointBoundary(const foamField &fieldName, label boundaryIndex) const;

    template<class type, label nComponents>
    [[nodiscard]] vtkSmartPointer<vtkFloatArray> convertEmptyVolBoundary(const foamField &fieldName, const fvPatch& boundaryPatch) const;

    template<class type>
    [[nodiscard]] vtkSmartPointer<vtkFloatArray> convertEmptyPointBoundary(const foamField &fieldName, const fvPatch& boundaryPatch) const;

    template<class T>
    [[nodiscard]] scalarMinMax calculateVectorMinMax(const foamField &fieldName, const T &myField) const;

    [[nodiscard]] scalarMinMax calculateDomainRangeForVectorField(const foamField &fieldName) const;
    [[nodiscard]] scalarMinMax calculateDomainRangeForPointVectorField(const foamField &fieldName) const;
    [[nodiscard]] scalarMinMax calculateDomainRangeForScalarField(const foamField &fieldName) const;
    [[nodiscard]] scalarMinMax calculateDomainRangeForPointScalarField(const foamField &fieldName) const;
    [[nodiscard]] static bool shouldFieldBeListed(const word &field);
};

} // End namespace functionObjects::runTimeVis
} // End namespace Foam
