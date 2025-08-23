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
    (c) 2020-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "pvdWriter.H"

#include "db/Time/Time.H"
#include "infos/runTimeVisualisationInfo.H"
#include "scene/items.H"
#include "types/compressionLevelType.H"

// VTK includes
#include "engysPVDWriter.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "Utils/basicProfiler.H"

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

PVDWriter::PVDWriter
    (
        const std::string& name,
        const RunTimeVisualisationInfo& rtppInfo,
        const SceneInfo& sceneInfo,
        vtkMultiProcessController* controller,
        Items &items
    )
    :
    items_(items)
{
    const PvdInfo& pvdInfo = rtppInfo.getPvdInfo();
    if (pvdInfo.exportPvd)
    {
        this->pvdWriter_ = vtkSmartPointer<engysPVDWriter>::New();
        this->pvdWriter_->SetController(controller);
        string fileName = rtppInfo.getOutputFolder() / name + ".pvd";
        this->pvdWriter_->SetFileName(fileName.c_str());
        this->pvdWriter_->SetDataModeToAppended();
        this->pvdWriter_->EncodeAppendedDataOff();
        this->pvdWriter_->SetWriteMultipleFilesInParallel(pvdInfo.oneFilePerProcess);
        this->pvdWriter_->SetRemoveOldObsoleteFiles(pvdInfo.removeOldObsoleteFiles);
        this->allowableFields_ = pvdInfo.getAllowedFieldsForPVD(sceneInfo);
        this->writePointData_ = pvdInfo.writePointData;

        switch (pvdInfo.compressionLevel.getValue())
        {
            case CompressionLevelType::FASTEST:
                pvdWriter_->SetCompressorTypeToNone();
                pvdWriter_->SetCompressionLevel(0);
                break;
            case CompressionLevelType::FAST:
                pvdWriter_->SetCompressorTypeToLZ4();
                pvdWriter_->SetCompressionLevel(3);
                break;
            default:
            case CompressionLevelType::BALANCED:
                pvdWriter_->SetCompressorTypeToLZ4();
                pvdWriter_->SetCompressionLevel(6);
                break;
            case CompressionLevelType::SMALL:
                pvdWriter_->SetCompressorTypeToZLib();
                pvdWriter_->SetCompressionLevel(6);
                break;
            case CompressionLevelType::SMALLEST:
                pvdWriter_->SetCompressorTypeToZLib();
                pvdWriter_->SetCompressionLevel(9);
                break;
        }
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void PVDWriter::exportIfRequired(const Time &currentTime)
{
    if (this->pvdWriter_)
    {
        vtkSmartPointer<vtkDataSet> itemData = this->items_.createDataSetWithAllData(
            currentTime.timeIndex(),
            currentTime.value());

        const scalar currentTimeValue = currentTime.time().value();
        const scalar userCurrentTimeValue = currentTime.time().timeToUserTime(currentTimeValue);
        const word currentTimeName = Time::timeName(userCurrentTimeValue);

        if (!writePointData_)
        {
            itemData->GetPointData()->Initialize();
        }
        if (!this->allowableFields_.hasAllFieldsMarker())
        {
            removeUnwantedFields(itemData->GetPointData(), allowableFields_.getPointFoamFieldsSet());
            removeUnwantedFields(itemData->GetCellData(), allowableFields_.getCellFoamFieldsSet());
        }

        this->pvdWriter_->SetCurrentTimeName(currentTimeName.c_str());
        this->pvdWriter_->SetInputData(itemData);
        basicProfiler::addMeasurePoint("Calling Write()");
        this->pvdWriter_->Write();
        this->pvdWriter_->RemoveAllInputs();
    }
}

void PVDWriter::removeUnwantedFields(vtkDataSetAttributes *data, std::set<std::string> wantedFields)
{
    std::vector<std::string> arraysToRemove;
    for (int i = 0; i < data->GetNumberOfArrays(); i++)
    {
        vtkDataArray* array = data->GetArray(i);
        std::string field(array->GetName());
        if (wantedFields.find(field) == wantedFields.end())
        {
            arraysToRemove.emplace_back(array->GetName());
        }
    }
    for (const std::string& arrayToRemove : arraysToRemove)
    {
        data->RemoveArray(arrayToRemove.c_str());
    }
}

} // End namespace Foam

// ************************************************************************* //
