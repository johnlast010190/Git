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
    (c) 2015-2019 OpenCFD Ltd.
    (c) 2022-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fieldSamplingDataSetProvider.H"

#include "dataStructs/objects/fieldSamplingObjectData.H"
#include "Utils/boundsUtils.H"

#include "engysResampleToImage.h"
#include "vtkImageData.h"
#include "vtkMultiProcessController.h"

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

FieldSamplingDataSetProvider::FieldSamplingDataSetProvider(
    const std::string &name,
    const FieldSamplingObjectData &dictData
)
    : ItemDataSetProvider(name)
{
    fieldSamplingFilter_ = vtkSmartPointer<engysResampleToImage>::New();
    scalar bounds[6];
    boundsUtils::scalarArrayFromBoundBox(bounds, dictData.samplingBounds);
    fieldSamplingFilter_->SetSamplingBounds(bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

    for (const word& field : dictData.desiredFields)
    {
        fieldSamplingFilter_->AddField(field.c_str());
    }

    fieldSamplingFilter_->SetSamplingDimensions(dictData.elementsX, dictData.elementsY, dictData.elementsZ);
    fieldSamplingFilter_->UseInputBoundsOff();
    fieldSamplingFilter_->SetController(vtkMultiProcessController::GetGlobalController());

    switch(dictData.toleranceType.getValue())
    {
        case runTimeVis::FieldSamplingToleranceType::Value::OFF:
            fieldSamplingFilter_->SetToleranceTypeToOff();
            break;
        case runTimeVis::FieldSamplingToleranceType::Value::ABSOLUTE:
            fieldSamplingFilter_->SetToleranceTypeToAbsolute();
            break;
        case runTimeVis::FieldSamplingToleranceType::Value::BOUNDS_RELATIVE:
            fieldSamplingFilter_->SetToleranceTypeToBoundsRelative();
            break;
        case runTimeVis::FieldSamplingToleranceType::Value::SPACING_RELATIVE:
            fieldSamplingFilter_->SetToleranceTypeToSpacingRelative();
            break;
    }
    fieldSamplingFilter_->SetTolerance(dictData.toleranceValue);
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void FieldSamplingDataSetProvider::update(scalar currentTime)
{
    vtkDataSet *sourceDataSet = sources_.at(0)->getDataSetOutput(1);
    fieldSamplingFilter_->SetInputDataObject(sourceDataSet);
    fieldSamplingFilter_->Update();

    output_ = fieldSamplingFilter_->GetOutput();
}

void FieldSamplingDataSetProvider::addExtraRequirementsDueToUpstreamRequirements(
    ItemRequirements &ownRequirements,
    const ItemRequirements &upstreamRequirements
) const
{
    for (const std::string& field : upstreamRequirements.getRequiredFields().getCellFoamFieldsSet())
    {
        ownRequirements.addToRequiredFields(foamField::withPointAssociation(field));
    }
}

} // End namespace Foam


// ************************************************************************* //
