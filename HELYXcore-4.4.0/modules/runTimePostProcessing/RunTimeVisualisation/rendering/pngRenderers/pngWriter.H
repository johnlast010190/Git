/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : Dev
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
    Foam::functionObjects::runTimeVis::PngWriter

Description
    base class for the classes that render the PNGs and write it to a file

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "primitives/strings/fileName/fileName.H"
#include "rendering/rendererWriter.H"

#include "vtkImageData.h"
#include "vtkSmartPointer.h"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declarations
class vtkRenderWindow;
class vtkImageData;

namespace Foam::functionObjects::runTimeVis
{

struct RenderInfo;

/*---------------------------------------------------------------------------*\
                       Class PngWriter Declaration
\*---------------------------------------------------------------------------*/

class PngWriter : public RendererWriter
{
public:

    explicit PngWriter(const RenderInfo& renderInfo) : renderInfo(renderInfo) {};

    ~PngWriter() override = default;

protected:
    vtkSmartPointer<vtkImageData> render(vtkRenderWindow* window) const;
    static void writeToFile(vtkImageData* imageData, fileName pathToOutputFileWithoutExtension);

private:
    const RenderInfo& renderInfo;
    vtkSmartPointer<vtkImageData> renderCroppedImage(vtkRenderWindow *window) const;
    vtkSmartPointer<vtkImageData> renderFullImage(vtkRenderWindow *window) const;
    static bool checkImageSize(vtkImageData *imageData);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
