/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
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
    (c) 2022-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "vectorWidget.H"

#include "db/IOobject/IOobject.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"

#include "dataStructs/widgets/vectorWidgetData.H"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "engysVectorWidget.h"

namespace Foam::functionObjects::runTimeVis
{

std::unique_ptr<VectorWidget> VectorWidget::createWidgetIfVisible(const VectorWidgetData& data)
{
    if (data.visible)
    {
        return std::unique_ptr<VectorWidget>(new VectorWidget(data));
    }
    else
    {
        return {nullptr};
    }
}


VectorWidget::VectorWidget(const VectorWidgetData& data)
{
    vectorWidget_ = vtkSmartPointer<engysVectorWidget>::New();

    dataSourceType_ = data.dataSource;

    vectorWidget_->SetShowVector(data.vector3dData.showVector3D);
    vectorWidget_->SetShowResultantVector(data.vector3dData.showResultant);
    vectorWidget_->SetShowComponentVectors(data.vector3dData.showComponents);
    vectorWidget_->SetVectorScreenSize(data.vector3dData.size);
    vectorWidget_->SetShowVectorTitle(data.vector3dData.showTitle);
    vectorWidget_->SetVectorTitle(data.dataSource.getTitle().c_str());
    vectorWidget_->SetVectorScreenCoordinates(data.vector3dData.coordinates.x(), data.vector3dData.coordinates.y());

    vectorWidget_->SetShowText(data.vectorTextData.showVectorTextRep);
    vectorWidget_->SetLabelFormat(data.vectorTextData.labelFormat.c_str());
    vectorWidget_->SetFontSize(data.vectorTextData.fontSize);
    vectorWidget_->SetShowTextBox(data.vectorTextData.showTextBox);
    vectorWidget_->SetShowTextTitle(data.vectorTextData.showTitle);
    vectorWidget_->SetTextTitle(data.dataSource.getTitle().c_str());
    vectorWidget_->SetTextCoordinates(data.vectorTextData.coordinates.x(), data.vectorTextData.coordinates.y());
}

void VectorWidget::addToRenderer(vtkRenderer* renderer)
{
    if (!vectorWidget_->GetInteractor())
    {
        vectorWidget_->SetInteractor(renderer->GetRenderWindow()->GetInteractor());
        vectorWidget_->Set2DRenderer(renderer);
    }
    vectorWidget_->On();
}

void VectorWidget::removeFromRenderer()
{
    vectorWidget_->Off();
}

void VectorWidget::update(const Time& runTime)
{
    const dictionary& runTimeInfoDict = runTime.lookupObject<IOdictionary>(vectorWidgetKeys::RUNTIME_INFO_DICT);
    const dictionary& resultsDict = runTimeInfoDict.subDict(vectorWidgetKeys::RESULTS_SUBDICT);
    const dictionary& solverDict = getSolverDict(resultsDict);
    const dictionary& vectorDict = solverDict.subDict(vectorWidgetKeys::VECTOR_SUBDICT_KEY);

    setReferenceVector(vectorDict);

    setVectorToShow(vectorDict);
}

const dictionary& VectorWidget::getSolverDict(const dictionary& resultsDict)
{
    if (resultsDict.isDict(vectorWidgetKeys::FLOWSOLVER_SUBDICT))
    {
        return resultsDict.subDict(vectorWidgetKeys::FLOWSOLVER_SUBDICT);
    }
    else if (resultsDict.isDict(vectorWidgetKeys::INTERFOAM_SUBDICT))
    {
        return resultsDict.subDict(vectorWidgetKeys::INTERFOAM_SUBDICT);
    }
    else
    {
        FatalErrorInFunction << "Could not find flowSolver nor interFoam as subdicts of results."
            << " Please report this error to the developers."
            << abort(FatalError);
        return resultsDict;
    }
}

void VectorWidget::setReferenceVector(const dictionary &vectorDict)
{
    if (referenceVectorSet_)
    {
        return;
    }

    auto maxBodyAccelerationComponents = vectorDict.lookup<Vector<scalar>>(vectorWidgetKeys::MAX_FRAME_ACCELERATION_COMPONENTS_KEY);
    auto minBodyAccelerationComponents = vectorDict.lookup<Vector<scalar>>(vectorWidgetKeys::MIN_FRAME_ACCELERATION_COMPONENTS_KEY);
    auto g = vectorDict.lookup<Vector<scalar>>(vectorWidgetKeys::GRAVITY_ACCELERATION_KEY);

    Vector<scalar> maxComponents = {0, 0, 0};
    Vector<scalar> minComponents = {0, 0, 0};

    if (dataSourceType_.hasBodyAcceleration())
    {
        maxComponents = maxBodyAccelerationComponents;
        minComponents = minBodyAccelerationComponents;
    }
    if (dataSourceType_.hasGravity())
    {
        maxComponents -= g;
        minComponents -= g;
    }
    vectorWidget_->SetReferenceVectors(maxComponents[0], maxComponents[1], maxComponents[2], minComponents[0], minComponents[1], minComponents[2]);

    referenceVectorSet_ = true;
}

void VectorWidget::setVectorToShow(const dictionary &vectorDict)
{
    auto frameAcceleration = vectorDict.lookup<Vector<scalar>>(vectorWidgetKeys::FRAME_ACCELERATION_KEY);
    auto gravity = vectorDict.lookup<Vector<scalar>>(vectorWidgetKeys::GRAVITY_ACCELERATION_KEY);

    switch(dataSourceType_.getValue())
    {
        case VectorDataType::BODY_ACCELERATION:
            vectorWidget_->SetVectorToShow(frameAcceleration[0], frameAcceleration[1], frameAcceleration[2]);
            break;
        case VectorDataType::BODY_ACCELERATION_AND_GRAVITY:
        {
            Vector<scalar> relative = frameAcceleration - gravity;
            vectorWidget_->SetVectorToShow(relative[0], relative[1], relative[2]);
            break;
        }
        case VectorDataType::GRAVITY:
            vectorWidget_->SetVectorToShow(gravity[0], gravity[1], gravity[2]);
            break;
        default:
            FatalErrorInFunction << "Unknown frame acceleration type!" << abort(FatalError);
    }
}

// ************************************************************************* //
} // End namespace Foam
