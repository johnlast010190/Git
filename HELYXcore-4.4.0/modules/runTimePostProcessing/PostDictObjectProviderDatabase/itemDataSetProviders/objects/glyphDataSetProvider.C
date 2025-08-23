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
    (c) 2024-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "glyphDataSetProvider.H"

#include "engysGlyphs.h"
#include "vtkPointData.h"

#include "vtkMultiProcessController.h"

namespace Foam::functionObjects::runTimeVis
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

GlyphDataSetProvider::GlyphDataSetProvider(const std::string& name, const GlyphObjectData& dictData)
: ItemDataSetProvider(name)
{
    vectorsFilter_ = vtkSmartPointer<engysGlyphs>::New();
    vectorsFilter_->SetDistributionType(dictData.glyphDistribution.type.getValue());
    vectorsFilter_->SetGlyphType(dictData.glyphType.getValue());
    vectorsFilter_->SetMaxGlyphLength(dictData.maxGlyphLength);
    vectorsFilter_->SetMaximumNumberOfPoints(dictData.glyphDistribution.maximumNumberOfPoints);
    vectorsFilter_->SetOrientationFieldName(dictData.orientationField.c_str());
    vectorsFilter_->SetPointsInterval(dictData.glyphDistribution.pointsInterval);
    vectorsFilter_->SetSeed(dictData.glyphDistribution.seed);
    vectorsFilter_->SetScalingField(
            dictData.enableScalingByField,
            dictData.scalingField.getFoamName().c_str(),
            dictData.scalingField.isMagnitude()
    );
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void GlyphDataSetProvider::update(scalar vtkNotUsed(currentTime))
{
    vtkSmartPointer<vtkDataSet> inputData = sources_.at(0)->getDataSetOutput();

    vectorsFilter_->SetInputData(inputData);
    vectorsFilter_->Update();
    output_ = vectorsFilter_->GetOutput();
}

} // End namespace


// ************************************************************************* //
