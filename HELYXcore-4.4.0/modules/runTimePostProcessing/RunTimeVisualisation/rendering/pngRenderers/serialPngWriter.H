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
    Foam::functionObjects::runTimeVis::SerialPngWriter

Description
    writes the png file when only the master is rendering

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "rendering/pngRenderers/pngWriter.H"
#include "Utils/ParallelUtils.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declarations

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                       Class SerialPngWriter Declaration
\*---------------------------------------------------------------------------*/

class SerialPngWriter : public PngWriter
{
public:
    explicit SerialPngWriter(const RenderInfo& renderInfo) : PngWriter(renderInfo) {};

    void renderAndWrite(
        vtkRenderWindow *window,
        const RendererExtraData& extraData
    ) const override
    {
        if (ParallelUtils::isMaster())
        {
            vtkSmartPointer<vtkImageData> imageData = render(window);
            writeToFile(imageData, extraData.pathToOutputFileWithoutExtension);
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
