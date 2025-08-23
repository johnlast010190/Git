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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "boundarySurfaceDataSetProvider.H"

#include "vtkPlaneSource.h"

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void BoundarySurfaceDataSetProvider::initialiseSurface()
{
    if (!ParallelUtils::isMaster() || data_.elements[0] == 0 || data_.elements[1] == 0)
    {
        output_ = vtkSmartPointer<vtkPolyData>::New();
        return;
    }

    vtkNew<vtkPlaneSource> planeSource;
    planeSource->SetOrigin(data_.points[0].x(), data_.points[0].y(), data_.points[0].z());
    planeSource->SetPoint1(data_.points[1].x(), data_.points[1].y(), data_.points[1].z());
    planeSource->SetPoint2(data_.points[3].x(), data_.points[3].y(), data_.points[3].z());
    planeSource->SetXResolution(data_.elements[0]);
    planeSource->SetYResolution(data_.elements[1]);
    planeSource->Update();

    output_ = planeSource->GetOutput();
}

} // End namespace


// ************************************************************************* //
