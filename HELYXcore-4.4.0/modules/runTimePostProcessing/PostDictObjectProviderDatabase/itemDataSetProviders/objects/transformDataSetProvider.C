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

#include "transformDataSetProvider.H"

#include "vtkTransform.h"
#include "vtkTransformFilter.h"
#include "vtkAppendFilter.h"
#include "vtkObjectFactory.h"


namespace Foam::functionObjects::runTimeVis
{

#define PI static_cast<scalar>(3.14159265358979)

// These classes are intended to me moved to the evtk in the future
class TransformWrapper : public vtkTransform {
    private:
    bool original_ = true;

    public:
    TransformWrapper() : vtkTransform(), original_(true) {
    }

    static TransformWrapper* New();

    explicit TransformWrapper(const TransformWrapper* wrapper) {
        copyFrom(wrapper);
    }

    void copyFrom(const TransformWrapper* wrapper) {
        DeepCopy(const_cast<TransformWrapper*>(wrapper));
        original_ = wrapper->isOriginal();
    }

    [[nodiscard]] bool isOriginal() const {
        return original_;
    }

    void setOriginal(bool original) {
        original_ = original;
    }

    void RotateWXYZAroundOrigin(Foam::point origin, Foam::vector axis, Foam::scalar degrees) {
        original_ = false;
        Translate(-origin[0], -origin[1], -origin[2]);
        RotateWXYZ(degrees, axis[0], axis[1], axis[2]);
        Translate(origin[0], origin[1], origin[2]);
    }

    void MirrorOriginAxisOffset(Foam::point origin, Foam::vector axis, Foam::scalar offsetMeters) {
        original_ = false;
        Translate(-origin[0], -origin[1], -origin[2]);

        scalar mag = Foam::mag(axis);
        scalar normX = axis[0] / mag;
        scalar normY = axis[1] / mag;
        scalar normZ = axis[2] / mag;

        if (normX * normX >= normY * normY && normX * normX >= normZ * normZ) {
            scalar degrees = std::acos(normX)*180/PI;
            RotateWXYZ(-degrees, 0, -normZ, normY);
            Scale(-1, 1, 1);
            RotateWXYZ(degrees, 0, -normZ, normY);
        } else if (normY * normY >= normZ * normZ) {
            scalar degrees = std::acos(normY)*180/PI;
            RotateWXYZ(-degrees, normZ, 0, -normX);
            Scale(1, -1, 1);
            RotateWXYZ(degrees, normZ, 0, -normX);
        } else {
            scalar degrees = std::acos(normZ)*180/PI;
            RotateWXYZ(-degrees, -normY, normX, 0);
            Scale(1, 1, -1);
            RotateWXYZ(degrees, -normY, normX, 0);
        }

        Translate(normX * offsetMeters, normY * offsetMeters, normZ * offsetMeters);

        Translate(origin[0], origin[1], origin[2]);
    }

    void TranslateVector(Foam::vector distance) {
        original_ = false;
        Translate(distance[0], distance[1], distance[2]);
    }

    void ScaleOrigin(Foam::vector ratio, Foam::point origin) {
        original_ = false;
        Translate(-origin[0], -origin[1], -origin[2]);
        Scale(ratio[0], ratio[1], ratio[2]);
        Translate(origin[0], origin[1], origin[2]);
    }
};

vtkStandardNewMacro(TransformWrapper)


// These classes are intended to me moved to the evtk in the future
class TransformationsGroup {
private:
    std::vector<vtkSmartPointer<TransformWrapper>> transforms_;

public:
    TransformationsGroup() {
        vtkSmartPointer<TransformWrapper> initialTransform = vtkSmartPointer<TransformWrapper>::New();
        initialTransform->PostMultiply();
        transforms_.push_back(initialTransform);
    }

    void applyTransformationByCode(const Foam::functionObjects::runTimeVis::TransformationType& code, Foam::point point, Foam::vector vector, Foam::scalar offset, Foam::label copyCount, bool includeSource) {
        std::vector<vtkSmartPointer<TransformWrapper>> tempTransformsList;
        for (TransformWrapper* transform : transforms_) {

            if (includeSource) {
                vtkSmartPointer<TransformWrapper> transformCopy = vtkSmartPointer<TransformWrapper>::New();
                transformCopy->copyFrom(transform);
                tempTransformsList.push_back(transformCopy);
            }

            switch (code.getValue()) {
                case Foam::functionObjects::runTimeVis::TransformationType::CYLINDRICAL_ARRAY:
                    for (label i = 0; i < copyCount; i++) {
                        vtkSmartPointer<TransformWrapper> transformCopy = vtkSmartPointer<TransformWrapper>::New();
                        transformCopy->copyFrom(transform);
                        transformCopy->RotateWXYZAroundOrigin(point, vector, offset * static_cast<scalar>(i + 1));
                        tempTransformsList.push_back(transformCopy);
                    }
                    break;

                case Foam::functionObjects::runTimeVis::TransformationType::MIRROR:
                    transform->MirrorOriginAxisOffset(point, vector, offset);
                    tempTransformsList.emplace_back(transform);
                    break;

                case Foam::functionObjects::runTimeVis::TransformationType::TRANSLATE:
                    transform->TranslateVector(vector);
                    tempTransformsList.emplace_back(transform);
                    break;

                case Foam::functionObjects::runTimeVis::TransformationType::ROTATE:
                    transform->RotateWXYZAroundOrigin(point, vector, offset);
                    tempTransformsList.emplace_back(transform);
                    break;

                case Foam::functionObjects::runTimeVis::TransformationType::LINEAR_ARRAY:
                {
                    auto translation = Foam::vector(vector);
                    translation.normalise();
                    translation *= offset;
                    vtkSmartPointer<TransformWrapper> previousTransform = transform;
                    for (label i = 0; i < copyCount; i++) {
                        vtkSmartPointer<TransformWrapper> transformCopy = vtkSmartPointer<TransformWrapper>::New();
                        transformCopy->copyFrom(previousTransform);
                        transformCopy->TranslateVector(translation);
                        tempTransformsList.push_back(transformCopy);
                        previousTransform = transformCopy;
                    }
                    break;
                }

                case Foam::functionObjects::runTimeVis::TransformationType::SCALE:
                    transform->ScaleOrigin(vector, point);
                    tempTransformsList.emplace_back(transform);
                    break;

                case Foam::functionObjects::runTimeVis::TransformationType::NONE:
                    tempTransformsList.emplace_back(transform);
                    break;
            }
        }
        transforms_ = tempTransformsList;
    }

    std::vector<vtkSmartPointer<vtkTransform>> getResultingTransformationsList() {
        std::vector<vtkSmartPointer<vtkTransform>> tempTransformsList;
        for (TransformWrapper* transform : transforms_) {
            tempTransformsList.emplace_back(transform);
        }
        return tempTransformsList;
    }

    std::vector<vtkSmartPointer<vtkTransform>> getResultingTransforms(const Foam::functionObjects::runTimeVis::TransformObjectData& transformInfo) {
        for (const Foam::functionObjects::runTimeVis::SingleTransformationObjectData& transformation : transformInfo.transformations) {
            applyTransformationByCode(transformation.type, transformation.transformPoint, transformation.transformVector, transformation.transformOffset, transformation.transformCount, transformation.includeSource);
        }
        return getResultingTransformationsList();
    }
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

TransformDataSetProvider::TransformDataSetProvider(const std::string& name, const TransformObjectData& dictData)
: ItemDataSetProvider(name)
{
    TransformationsGroup group;
    std::vector<vtkSmartPointer<vtkTransform>> transformsVTK = group.getResultingTransforms(dictData);
    appendFilter_ = vtkSmartPointer<vtkAppendFilter>::New();

    for (auto & i : transformsVTK) {
        // We must recreate the filter each time, because the output it gives is linked to the filter,
        // which means that when the filter updates, all previous outputs it returned are updated as well
        vtkNew<vtkTransformFilter> transformFilter;
        transformFilters_.emplace_back(transformFilter);
        transformFilter->SetTransform(i);
        transformFilter->TransformAllInputVectorsOn();
        appendFilter_->AddInputConnection(transformFilter->GetOutputPort());
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void TransformDataSetProvider::update(scalar vtkNotUsed(currentTime))
{
    for (vtkTransformFilter* filter : transformFilters_)
    {
        filter->SetInputData(sources_.at(0)->getDataSetOutput());
        filter->Update();
    }
    appendFilter_->Update();
    output_ = appendFilter_->GetOutput();
}

} // End namespace


// ************************************************************************* //
