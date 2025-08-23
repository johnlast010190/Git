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

#include "cuttingSurfaceProviders.H"

#include "vtkBox.h"
#include "vtkTransform.h"

namespace Foam::functionObjects::runTimeVis
{

CuttingBoxProvider::CuttingBoxProvider(const CuttingBoxData& data)
{
    vtkSmartPointer<vtkBox> box = vtkSmartPointer<vtkBox>::New();

    box->SetBounds
    (
        data.position.x() - data.scale.x()/2.0,  // Min x
        data.position.x() + data.scale.x()/2.0,  // Max x
        data.position.y() - data.scale.y()/2.0,  // Min y
        data.position.y() + data.scale.y()/2.0,  // Max y
        data.position.z() - data.scale.z()/2.0,  // Min z
        data.position.z() + data.scale.z()/2.0   // Max z
    );

    if
    (
        data.rotation.x() != 0 ||
        data.rotation.y() != 0 ||
        data.rotation.z() != 0
    )
    {
        vtkSmartPointer<vtkTransform> t =
            vtkSmartPointer<vtkTransform>::New();
        // Note from the GUI:
        // This transformation need to have all the signs opposite respect
        // the others. Like the axes are inverted.
        t->Translate
        (
            data.position.x(),
            data.position.y(),
            data.position.z()
        );
        t->RotateY(-data.rotation.y());
        t->RotateX(-data.rotation.x());
        t->RotateZ(-data.rotation.z());
        t->Translate
        (
            -data.position.x(),
            -data.position.y(),
            -data.position.z()
        );
        box->SetTransform(t);
    }
    implicitFunction_ = box;
}

} // End namespace


// ************************************************************************* //
