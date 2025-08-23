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
    (c) 2020-2023 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "surfaceTransformationsData.H"

#include "types/surfaceTransformationType.H"
#include "postDict/postDictKeys.H"
#include "dataStructs/geometry/surfaceTransformations/surfaceScaleData.H"
#include "dataStructs/geometry/surfaceTransformations/surfaceTranslateData.H"
#include "dataStructs/geometry/surfaceTransformations/surfaceRotateData.H"

#include "vtkTransform.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations

namespace Foam::functionObjects::runTimeVis
{

SurfaceTransformationsData::SurfaceTransformationsData(const dictionary& geometryDict)
{
    readDict(geometryDict);
}

SurfaceTransformationsData::SurfaceTransformationsData(const SurfaceTransformationsData &t)
{
    transformations.clear();
    for (const auto& transformation : t.transformations)
    {
        transformations.emplace_back(transformation->copy());
    }
}

void SurfaceTransformationsData::readDict(const dictionary &geometryDict)
{
    transformations.clear();
    if (!geometryDict.found(surfaceTransformationKeys::TRANSFORMS_KEY))
    { return; }

    auto transformationDictList = geometryDict.lookup<List<dictionary>>(surfaceTransformationKeys::TRANSFORMS_KEY);
    for (const dictionary &transformationDict: transformationDictList)
    {
        auto type = transformationDict.lookup<SurfaceTransformationType>(
            surfaceTransformationKeys::TYPE_KEY
        );
        switch (type.getValue())
        {
            case SurfaceTransformationType::SCALE:
                transformations.emplace_back(
                    std::unique_ptr<SurfaceTransformationData>(
                        new SurfaceScaleData(
                            transformationDict
                        )));
                break;
            case SurfaceTransformationType::TRANSLATE:
                transformations.emplace_back(
                    std::unique_ptr<SurfaceTransformationData>(
                        new SurfaceTranslateData(
                            transformationDict
                        )));
                break;
            case SurfaceTransformationType::ROTATE:
                transformations.emplace_back(
                    std::unique_ptr<SurfaceTransformationData>(
                        new SurfaceRotateData(
                            transformationDict
                        )));
                break;
            case SurfaceTransformationType::UNKNOWN:
                break;
        }
    }
}

bool SurfaceTransformationsData::operator==(const SurfaceTransformationsData &other) const
{
    if (transformations.size() != other.transformations.size())
    { return false; }

    for (size_t i = 0; i < transformations.size(); i++)
    {
        if (transformations[i] != other.transformations[i])
        { return false; }
    }
    return true;
}

void SurfaceTransformationsData::computeAndAddToHash(size_t &hash) const
{
    for (const auto &transformation: transformations)
    {
        transformation->computeAndAddToHash(hash);
    }
}

vtkSmartPointer<vtkTransform> SurfaceTransformationsData::getVtkTransform() const
{
    vtkNew<vtkTransform> transform;
    transform->PostMultiply();
    for (const auto &transformation: transformations)
    {
        transformation->addToVtkTransform(transform);
    }
    return transform;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
