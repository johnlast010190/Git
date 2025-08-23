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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "logoWidget.H"

#include "dataStructs/widgets/logoWidgetData.H"

#include "global/etcFiles/etcFiles.H"
#include "include/OSspecific.H"

#include "vtkPNGReader.h"
#include "vtkActor2D.h"
#include "vtkImageMapper.h"
#include "vtkImageData.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"

static bool alreadyDisplayedWarning = false;

namespace Foam::functionObjects::runTimeVis
{

std::unique_ptr<LogoWidget> LogoWidget::createWidgetIfVisible(const LogoWidgetData& data)
{
    if (data.visible)
    {
        return std::unique_ptr<LogoWidget>(new LogoWidget());
    }
    else
    {
        return nullptr;
    }
}

void LogoWidget::initialise(int imageWidth)
{
    initialised_ = true;
    fileName imgName = findEtcFile("../../../GUI/img/Logo.png");

    if (word::null == imgName) {
        imgName = getEnv("HELYX_GUI_PATH")/"img"/"Logo.png";
        if (!isFile(imgName))
        {
            imgName = word::null;
        }
    }

    if (word::null == imgName)
    {
        if (!alreadyDisplayedWarning)
        {
            alreadyDisplayedWarning = true;
            Warning
                << "The Engys logo was not found, and will not be present "
                << "in generated images (this is normally the case with "
                << "self-compiled versions of HELYXcore)"
                << endl;
        }
    }
    else
    {
        vtkNew<vtkPNGReader> reader;
        reader->SetFileName(imgName.c_str());
        reader->Update();
        vtkSmartPointer<vtkImageData> imageData = reader->GetOutput();

        auto logoWidth = static_cast<label>(imageData->GetExtent()[1]);
        auto logoHeight = static_cast<label>(imageData->GetExtent()[3]);
        label bottomLeftCornerX = imageWidth - logoWidth - 40;
        label bottomLeftCornerY = 40;
        label topRightCornerX = logoWidth + 1;
        label topRightCornerY = logoHeight + 1;

        vtkNew<vtkImageMapper> mapper;
        mapper->SetInputData(imageData);
        mapper->SetColorWindow(255);
        mapper->SetColorLevel(127.5);

        actor_ = vtkSmartPointer<vtkActor2D>::New();
        actor_->SetPosition(bottomLeftCornerX, bottomLeftCornerY);
        actor_->SetPosition2(topRightCornerX, topRightCornerY);
        actor_->SetMapper(mapper);
    }
}

void LogoWidget::addToRenderer(vtkRenderer* renderer)
{
    if (!initialised_)
    {
        initialise(renderer->GetRenderWindow()->GetSize()[0]);
    }
    if (actor_)
    {
        renderer->AddActor2D(actor_);
    }
}

void LogoWidget::removeFromRenderer(vtkRenderer* renderer)
{
    if (actor_)
    {
        renderer->RemoveActor2D(actor_);
    }
}

// ************************************************************************* //
} // End namespace
