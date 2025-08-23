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

Class
    Foam::functionObjects::runTimeVis::SerialX3dWriter

Description
    writes the x3d file with only the master

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "primitives/strings/fileName/fileName.H"
#include "include/OSspecific.H"

#include "Utils/vtkDataHandlingTools.H"
#include "Utils/ParallelUtils.H"

#include "vtkSmartPointer.h"
#include "vtkRenderWindow.h"
#include "vtkX3DExporter.h"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declarations

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                       Class SerialX3dWriter Declaration
\*---------------------------------------------------------------------------*/

class SerialX3dWriter
{
public:
    static void writeX3D(vtkRenderWindow* window, vtkRenderer* renderer, fileName pathToOutputFileWithoutExtension)
    {
        if (ParallelUtils::isMaster())
        {
            vtkSmartPointer<vtkX3DExporter> exporter =
                vtkSmartPointer<vtkX3DExporter>::New();
            exporter->SetInput(window);
            exporter->SetActiveRenderer(renderer);

            fileName x3dOutputFile = pathToOutputFileWithoutExtension.ext("x3d");
            mkDir(x3dOutputFile.path());
            exporter->SetFileName(x3dOutputFile.c_str());

            Info<< "Generating file " << x3dOutputFile << endl;

            exporter->Update();
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
