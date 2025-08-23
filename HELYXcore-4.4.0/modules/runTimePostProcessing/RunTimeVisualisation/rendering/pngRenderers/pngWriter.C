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
    (c) 2020-2021 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::PngWriter

Description
    base class for the classes that render the PNGs and write it to a file

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#include "pngWriter.H"

#include "infos/renderInfo.H"

#include "Utils/vtkDataHandlingTools.H"
#include "Utils/ParallelUtils.H"

#include "vtkSmartPointer.h"
#include "vtkRenderWindow.h"
#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"
#include "vtkRendererCollection.h"
#include "engysWindowToCroppedImage.h"

#include "include/OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declarations

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                       Class PngWriter Declaration
\*---------------------------------------------------------------------------*/

vtkSmartPointer<vtkImageData> PngWriter::render(vtkRenderWindow* window) const
{
    if (renderInfo.cropToContents)
    {
        return renderCroppedImage(window);
    }
    else
    {
        return renderFullImage(window);
    }
}

vtkSmartPointer<vtkImageData> PngWriter::renderCroppedImage(vtkRenderWindow *window) const
{
    vtkNew<engysWindowToCroppedImage> windowToImageFilter;
    windowToImageFilter->SetCenterImageParts(renderInfo.centerToCroppedContents);
    windowToImageFilter->SetRenderWindow(window);
    vtkRendererCollection* renderers = window->GetRenderers();
    renderers->InitTraversal();
    while (vtkRenderer* renderer = renderers->GetNextItem())
    {
        if (renderer->GetLayer() == 0)
        {
            windowToImageFilter->SetMainRenderer(renderer);
            break;
        }
    }
    windowToImageFilter->Update();
    return windowToImageFilter->GetOutput();
}

vtkSmartPointer<vtkImageData> PngWriter::renderFullImage(vtkRenderWindow *window) const
{
    window->Render();
    vtkNew<vtkWindowToImageFilter> windowToImageFilter;
    windowToImageFilter->SetInput(window);
    windowToImageFilter->ShouldRerenderOff();
    if (renderInfo.transparentBackground)
    {
        windowToImageFilter->SetInputBufferTypeToRGBA();
    }
    else
    {
        windowToImageFilter->SetInputBufferTypeToRGB();
    }
    windowToImageFilter->Update();
    return windowToImageFilter->GetOutput();
}

void PngWriter::writeToFile(vtkImageData* imageData, fileName pathToOutputFileWithoutExtension)
{
    if (ParallelUtils::isMaster() && checkImageSize(imageData))
    {
        vtkNew<vtkPNGWriter> writer;

        fileName pngOutputFile = pathToOutputFileWithoutExtension.ext("png");
        writer->SetFileName(pngOutputFile.c_str());
        writer->SetInputData(imageData);

        mkDir(pngOutputFile.path());

        Info<< "Generating image " << pngOutputFile << endl;
        writer->Write();
    }
}

bool PngWriter::checkImageSize(vtkImageData *imageData)
{
    int extent[6];
    imageData->GetExtent(extent);
    int width = extent[1] - extent[0] + 1;
    int height = extent[3] - extent[2] + 1;
    if (width <= 0 || height <= 0)
    {
        Warning << "Trying to write an image with invalid size (" << width << " x " << height << "). "
        << "Check if the scene has been defined correctly. The image will not be written." << endl;
        return false;
    }
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
