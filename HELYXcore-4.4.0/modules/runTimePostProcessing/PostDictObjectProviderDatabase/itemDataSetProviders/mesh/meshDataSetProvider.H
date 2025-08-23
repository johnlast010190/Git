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
    (c) 2022-2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::MeshDataSetProvider

Description
    Base class for the mesh providers, which computes their vtk dataset

SourceFiles
    meshDataSetProvider.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "itemDataSetProviders/itemDataSetProvider.H"
#include "baseClasses/meshAndFields.H"

#include "primitives/ranges/MinMax/MinMax.H"
#include "volMesh/volMesh.H"
#include "vtk/adaptor/foamVtkVtuAdaptor.H"

#include <unordered_map>
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


class vtkPointDataToCellData;
class vtkCellDataToPointData;

// Forward Declarations
namespace Foam
{
class fvMesh;
}

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                   Class MeshDataSetProvider Declaration
\*---------------------------------------------------------------------------*/

class MeshDataSetProvider : public ItemDataSetProvider
{
protected:

    // Protected data
    const string region_;

    const MeshAndFields meshAndFields_;

    vtk::vtuAdaptor adaptor_;

public:
    MeshDataSetProvider(const std::string &name, string region, const MeshAndFields& meshAndFields);
    ~MeshDataSetProvider() override = 0;

    [[nodiscard]] string getName() override
    { return name_; };

protected:
    vtkSmartPointer<vtkDataArray> getVolmeshFieldDataArray(const foamField &fieldName, bool& isPointArray);
    vtkSmartPointer<vtkDataSet> setGhostCells(vtkDataSet *mesh, bool isInternalMesh = false);
    void addInterpolatedCellAndPointData(vtkDataSet *dataset);

    [[nodiscard]] const FoamFields &getRequiredFields();

    vtkSmartPointer<vtkPointData> interpolateCellToPoint(vtkDataSet *dataset);
    static vtkSmartPointer<vtkCellData> interpolatePointToCell(vtkDataSet *dataset);
    static void addInterpolatedData(vtkDataSetAttributes *baseData, vtkDataSetAttributes *interpolatedData);

private:

    FoamFields updatedRequiredFields_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
