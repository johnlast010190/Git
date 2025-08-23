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

#include "thresholdDataSetProvider.H"

#include "engysThreshold.h"

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ThresholdDataSetProvider::ThresholdDataSetProvider(const std::string& name, const ThresholdObjectData& dictData)
        : ItemDataSetProvider(name)
{
    thresholdFilter_ = vtkSmartPointer<engysThreshold>::New();

    thresholdFilter_->SetIsMagnitude(dictData.field.isMagnitude());
    thresholdFilter_->SetIsPointAssociation(dictData.field.isPointAssociation());
    thresholdFilter_->SetNameOfFieldToProcess(dictData.field.getFoamName().c_str());

    thresholdFilter_->SetAllPointsCriterion(dictData.allPointsCriterion);
    thresholdFilter_->SetMinThreshold(dictData.minThreshold);
    thresholdFilter_->SetMaxThreshold(dictData.maxThreshold);

    // Following comment from the GUI:
    // NOTICE: "invert" not present in our vtk code, to be added when it's time
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void ThresholdDataSetProvider::update(scalar vtkNotUsed(currentTime))
{
    vtkDataSet* sourceDataSet = sources_.at(0)->getDataSetOutput();
    thresholdFilter_->SetInputData(sourceDataSet);
    thresholdFilter_->Update();
    output_ = thresholdFilter_->GetOutput();
}

} // End namespace

// ************************************************************************* //
