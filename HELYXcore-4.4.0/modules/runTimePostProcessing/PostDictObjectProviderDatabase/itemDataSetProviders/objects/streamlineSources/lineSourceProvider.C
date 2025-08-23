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

#include "streamlineSourceProviders.H"

#include "vtkLineSource.h"

namespace Foam::functionObjects::runTimeVis
{

LineSourceProvider::LineSourceProvider(const LineSourceData& data)
{
    vtkNew<vtkLineSource> seeds;
    seeds->SetPoint1
    (
        data.point1.x(),
        data.point1.y(),
        data.point1.z()
    );
    seeds->SetPoint2
    (
        data.point2.x(),
        data.point2.y(),
        data.point2.z()
    );
    seeds->SetResolution(data.numberOfPoints);
    seeds->Update();

    seedPoints_ = seeds->GetOutput();
}

vtkSmartPointer<vtkPolyData> LineSourceProvider::getSeedPoints()
{
    return seedPoints_;
}

} // End namespace


// ************************************************************************* //
