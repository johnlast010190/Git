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
    Foam::functionObjects::runTimeVis::EdfWriter

\*---------------------------------------------------------------------------*/

#include "edfWriter.H"

#include "infos/renderInfo.H"

#include "Utils/vtkDataHandlingTools.H"
#include "Utils/ParallelUtils.H"

#include "rendering/pngRenderers/parallelPngWriter.H"

#include "vtkSmartPointer.h"
#include "vtkRenderWindow.h"
#include "engysEDFOpenGLRenderer.h"
#include "engysEDF.h"
#include "engysEDFWriter.h"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declarations

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                       Class EdfWriter Declaration
\*---------------------------------------------------------------------------*/

EdfWriter::EdfWriter(const RenderInfo &renderInfo, engysEDFOpenGLRenderer* renderer)
    : renderer(renderer), edfCompressionLevel(renderInfo.edfCompressionLevel)
{
    if (renderInfo.exportPng)
    {
        this->delegateWriter = std::make_unique<ParallelPngWriter>(renderInfo);
    }
}

void EdfWriter::renderAndWrite(
    vtkRenderWindow *window,
    const RendererExtraData& extraData
) const
{
    if (this->delegateWriter)
    {
        this->delegateWriter->renderAndWrite(window, extraData);
    }
    else
    {
        window->Render();
    }
    vtkSmartPointer<engysEDF> edf = this->renderer->GetEDFResult(
        extraData.currentTime,
        extraData.sceneName.c_str(),
        extraData.cameraName.c_str());
    writeToFile(edf, extraData.pathToOutputFileWithoutExtension);
}

void EdfWriter::writeToFile(engysEDF* edf, fileName pathToOutputFileWithoutExtension) const
{
    if (ParallelUtils::isMaster() && edf->IsEDFValid())
    {
        vtkNew<engysEDFWriter> writer;

        fileName pngOutputFile = pathToOutputFileWithoutExtension.ext("edf");
        writer->SetFileName(pngOutputFile.c_str());
        writer->SetDataModeToAppended();
        writer->EncodeAppendedDataOff();
        SetupCompressionLevel(writer);
        writer->SetInputData(edf);

        mkDir(pngOutputFile.path());

        Info << "Generating edf " << pngOutputFile << endl;
        writer->Write();
    }
}

void EdfWriter::SetupCompressionLevel(engysEDFWriter *writer) const
{
    switch (edfCompressionLevel)
    {
        case FASTEST:
            writer->PNGCompressionOff();
            writer->SetCompressorTypeToLZ4();
            writer->SetCompressionLevel(1);
            break;
        case FAST:
            writer->PNGCompressionOff();
            writer->SetCompressorTypeToZLib();
            writer->SetCompressionLevel(3);
            break;
        case BALANCED:
            writer->PNGCompressionOff();
            writer->SetCompressorTypeToLZ4();
            writer->SetCompressionLevel(4);
            break;
        case SMALL:
            writer->PNGCompressionOn();
            writer->SetCompressionLevel(1);
            break;
        case SMALLEST:
            writer->PNGCompressionOn();
            writer->SetCompressionLevel(6);
            break;
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
