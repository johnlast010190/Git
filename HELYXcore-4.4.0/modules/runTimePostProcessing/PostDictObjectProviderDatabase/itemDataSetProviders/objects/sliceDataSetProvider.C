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

#include "sliceDataSetProvider.H"

#include "storage/referenceFrames.H"

#include "engysSliceDataSet.h"
#include "Utils/Utils.H"
#include "vtkTransform.h"

namespace Foam::functionObjects::runTimeVis
{

SliceDataSetProvider::SliceDataSetProvider
(
        const std::string& name,
        const SliceObjectData& dictData,
        const std::shared_ptr<const CuttingSurfaceProvider>& cutProvider,
        const ReferenceFrames& referenceFrames
)
:
        ItemDataSetProvider(name),
        cutProvider_(cutProvider)
{
    implicitFunction_ = cutProvider->getImplicitFunction();

    sliceAlgorithm_ = vtkSmartPointer<engysSliceDataSet>::New();
    sliceAlgorithm_->SetCrinkle(dictData.crinkle);
    if (!dictData.offsets.empty())
    {
        sliceAlgorithm_->SetNumberOfContours(dictData.offsets.size());
        int i = 0;
        for (const scalar& offset : dictData.offsets)
        {
            sliceAlgorithm_->SetValue(i++, offset);
        }
    }


    if (!dictData.referenceFrame.empty())
    {
        this->referenceFrame_ = &referenceFrames.getReferenceFrame(dictData.referenceFrame);
    }
    else
    {
        this->referenceFrame_ = nullptr;
        sliceAlgorithm_->SetSliceFunction(implicitFunction_);
    }
}

void SliceDataSetProvider::update(scalar vtkNotUsed(currentTime))
{
    vtkSmartPointer<vtkTransform> originalTransform;
    if (referenceFrame_)
    {
        originalTransform = vtkTransform::SafeDownCast(implicitFunction_->GetTransform());

        vtkNew<vtkTransform> modifiedTransform;
        if (originalTransform)
        {
            modifiedTransform->DeepCopy(originalTransform);
        }
        else
        {
            modifiedTransform->Identity();
        }
        modifiedTransform->PreMultiply();
        vtkSmartPointer<vtkMatrix4x4> referenceFrameMatrix = referenceFrame_->toTransformMatrix();
        referenceFrameMatrix->Invert();
        modifiedTransform->Concatenate(referenceFrameMatrix);
        implicitFunction_->SetTransform(modifiedTransform);
        sliceAlgorithm_->SetSliceFunction(implicitFunction_);
    }

    vtkDataSet* sourceDataSet = sources_.at(0)->getDataSetOutput();
    sliceAlgorithm_->SetInputDataObject(sourceDataSet);
    sliceAlgorithm_->Update();
    output_ = vtkDataSet::SafeDownCast(sliceAlgorithm_->GetOutputDataObject(0));

    if (referenceFrame_)
    {
        implicitFunction_->SetTransform(originalTransform);
    }
}

} // End namespace

// ************************************************************************* //
