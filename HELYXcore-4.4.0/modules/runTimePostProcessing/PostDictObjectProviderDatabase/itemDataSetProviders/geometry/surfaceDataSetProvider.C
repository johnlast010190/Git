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

\*---------------------------------------------------------------------------*/

#include "surfaceDataSetProvider.H"

#include "engysOutlineFilter.h"
#include "vtkBoundingBox.h"
#include "vtkCubeSource.h"
#include "vtkMultiProcessController.h"
#include "vtkTransformFilter.h"
#include "vtkTransform.h"

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

SurfaceDataSetProvider::SurfaceDataSetProvider(const std::string &name, const SurfaceTransformationsData& transformationsData)
: ItemDataSetProvider(name), transformationsData_(transformationsData), originalBounds_()
{
    vtkBoundingBox box;
    box.GetBounds(originalBounds_);
    outline_ = nullptr;
}

SurfaceDataSetProvider::SurfaceDataSetProvider(const std::string &name)
    : ItemDataSetProvider(name), transformationsData_(), originalBounds_()
{
    vtkBoundingBox box;
    box.GetBounds(originalBounds_);
    outline_ = nullptr;
}

SurfaceDataSetProvider::~SurfaceDataSetProvider() = default;

void SurfaceDataSetProvider::update(scalar currentTime)
{
    if (!initialised_)
    {
        initialiseSurface();

        vtkSmartPointer<vtkDataSet> surface = getDataSetOutput();
        if (surface->GetNumberOfPoints() > 0)
        {
            surface->GetBounds(originalBounds_);
        }

        transformSurface();

        initialised_ = true;
    }
}

void SurfaceDataSetProvider::transformSurface()
{
    if (!transformationsData_.transformations.empty())
    {
        vtkNew<vtkTransformFilter> transformFilter;
        transformFilter->SetInputData(output_);
        transformFilter->SetTransform(transformationsData_.getVtkTransform());
        transformFilter->Update();
        output_ = transformFilter->GetOutput();
    }
}

vtkSmartPointer<vtkPolyData> SurfaceDataSetProvider::getDataSetOutline()
{
    if (!outline_)
    {
        vtkSmartPointer<vtkPolyData> polyData;

        vtkBoundingBox box(originalBounds_);
        if (box.IsValid())
        {
            vtkNew<vtkCubeSource> cubeSource;
            cubeSource->SetBounds(originalBounds_);
            cubeSource->Update();
            polyData = cubeSource->GetOutput();
        }
        else
        {
            polyData = vtkSmartPointer<vtkPolyData>::New();
        }

        vtkNew<engysOutlineFilter> outlineFilter;
        outlineFilter->AllReduceOff();
        outlineFilter->SetController(vtkMultiProcessController::GetGlobalController());

        outlineFilter->SetInputData(polyData);
        outlineFilter->Update();
        outline_ = outlineFilter->GetOutput();

        if (!transformationsData_.transformations.empty())
        {
            vtkNew<vtkTransformFilter> transformFilter;
            transformFilter->SetInputData(outline_);
            transformFilter->SetTransform(transformationsData_.getVtkTransform());
            transformFilter->Update();
            outline_ = transformFilter->GetPolyDataOutput();
        }
    }

    return outline_;

}


} // End namespace Foam


// ************************************************************************* //
