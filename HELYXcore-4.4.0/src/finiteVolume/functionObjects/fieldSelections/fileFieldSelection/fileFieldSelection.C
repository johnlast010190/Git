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
    (c) 2017-19 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "functionObjects/fieldSelections/fileFieldSelection/fileFieldSelection.H"
#include "db/objectRegistry/objectRegistry.H"
#include "volMesh/volMesh.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchField.H"
#include "surfaceMesh/surfaceMesh.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchField.H"
#include "meshes/pointMesh/pointMesh.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchField.H"
#include "fields/UniformDimensionedFields/UniformDimensionedField.H"

void Foam::functionObjects::fileFieldSelection::addInternalFieldTypes
(
    DynamicList<fieldInfo>& set
) const
{
    const fvMesh& mesh = static_cast<const fvMesh&>(obr_);

    const IOobjectList allObjects(mesh, mesh.time().timeName());

    addFromFile<DimensionedField<scalar, volMesh>>(allObjects, set);
    addFromFile<DimensionedField<vector, volMesh>>(allObjects, set);
    addFromFile<DimensionedField<sphericalTensor, volMesh>>(allObjects, set);
    addFromFile<DimensionedField<symmTensor, volMesh>>(allObjects, set);
    addFromFile<DimensionedField<tensor, volMesh>>(allObjects, set);
}


void Foam::functionObjects::fileFieldSelection::addUniformFieldTypes
(
    DynamicList<fieldInfo>& set
) const
{
    const fvMesh& mesh = static_cast<const fvMesh&>(obr_);

    const IOobjectList allObjects(mesh, mesh.time().timeName());

    addFromFile<UniformDimensionedField<scalar>>(allObjects, set);
    addFromFile<UniformDimensionedField<vector>>(allObjects, set);
    addFromFile<UniformDimensionedField<sphericalTensor>>(allObjects, set);
    addFromFile<UniformDimensionedField<symmTensor>>(allObjects, set);
    addFromFile<UniformDimensionedField<tensor>>(allObjects, set);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fileFieldSelection::fileFieldSelection
(
    const objectRegistry& obr,
    const bool includeComponents
)
:
    fieldSelection(obr, includeComponents)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fileFieldSelection::updateSelection()
{
    List<fieldInfo> oldSet(std::move(selection_));

    DynamicList<fieldInfo> newSelection(oldSet.size());

    // Geometric fields
    addGeoFieldTypes<fvPatchField, volMesh>(newSelection);
    addGeoFieldTypes<fvsPatchField, surfaceMesh>(newSelection);
    addGeoFieldTypes<pointPatchField, pointMesh>(newSelection);

    // Internal fields
    addInternalFieldTypes(newSelection);

    // Uniform fields
    addUniformFieldTypes(newSelection);

    selection_.transfer(newSelection);

    (void)fieldSelection::checkSelection();

    return selection_ != oldSet;
}


// ************************************************************************* //
