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
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "surfacePlaneDataSetProvider.H"

#include "vtkPlaneSource.h"
#include "vtkTriangleFilter.h"

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void SurfacePlaneDataSetProvider::initialiseSurface()
{
    if (!ParallelUtils::isMaster())
    {
        output_ = vtkSmartPointer<vtkPolyData>::New();
        return;
    }

    vtkNew<vtkPlaneSource> planeSource;
    planeSource->SetOrigin(0, 0, 0);

    scalar diagonal = getPlaneDiagonal();
    planeSource->SetPoint1(diagonal, 0, 0);
    planeSource->SetPoint2(0, diagonal, 0);
    planeSource->SetCenter(data_.basePoint.x(), data_.basePoint.y(), data_.basePoint.z());
    planeSource->SetNormal
        (
            data_.normalVector.x(),
            data_.normalVector.y(),
            data_.normalVector.z()
        );

    vtkNew<vtkTriangleFilter> filter;
    filter->SetInputConnection(planeSource->GetOutputPort());
    filter->Update();
    output_ = filter->GetOutput();
}

scalar SurfacePlaneDataSetProvider::getPlaneDiagonal() const
{
    if (std::isnan(data_.diagonal))
    {
        return static_cast<scalar>(2.0) * data_.meshBounds.mag();
    }
    else
    {
        return data_.diagonal;
    }
}

} // End namespace Foam


// ************************************************************************* //
