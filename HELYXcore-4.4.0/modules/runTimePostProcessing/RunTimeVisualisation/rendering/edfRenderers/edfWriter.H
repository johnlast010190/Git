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

Description
    Renders the Engys Dynamic Fields EDFs and writes them to files

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "primitives/strings/fileName/fileName.H"
#include "rendering/rendererWriter.H"
#include "infos/renderInfo.H"

#include "vtkImageData.h"
#include "vtkSmartPointer.h"

#include "include/OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declarations
class vtkRenderWindow;
class vtkImageData;
class engysEDFOpenGLRenderer;
class engysEDF;
class engysEDFWriter;

namespace Foam::functionObjects::runTimeVis
{

class PngWriter;

/*---------------------------------------------------------------------------*\
                       Class EdfWriter Declaration
\*---------------------------------------------------------------------------*/

class EdfWriter : public RendererWriter
{
public:

    EdfWriter(const RenderInfo& renderInfo, engysEDFOpenGLRenderer* renderer);

    void renderAndWrite(
        vtkRenderWindow *window,
        const RendererExtraData& extraData
    ) const override;

    ~EdfWriter() override = default;

protected:
    void writeToFile(engysEDF* edf, fileName pathToOutputFileWithoutExtension) const;

private:
    std::unique_ptr<PngWriter> delegateWriter;
    vtkSmartPointer<engysEDFOpenGLRenderer> renderer;
    EDFCompressionLevel edfCompressionLevel;
    void SetupCompressionLevel(engysEDFWriter *writer) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
