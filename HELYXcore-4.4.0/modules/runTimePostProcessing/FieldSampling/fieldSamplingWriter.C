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

#include "fieldSamplingWriter.H"

#include "fieldSamplingData.H"
#include "compileOptions.H"

#include "Utils/boundsUtils.H"
#include "postDictObjectProviderDatabase.H"

#include "engysPVDWriter.h"
#include "engysResampleToImage.h"
#include "vtkMultiProcessController.h"
#include "vtkErrorCode.h"

#if BASIC_PROFILE == 1
#include "Utils/basicProfiler.H"
#define ADD_MEASURE_POINT(label) runTimeVis::basicProfiler::addMeasurePoint(label);
#else
#define ADD_MEASURE_POINT(label)
#endif


namespace Foam::functionObjects::fieldSample
{

FieldSamplingWriter::FieldSamplingWriter(
    const FieldSamplingData &data,
    const Switch &log
)
    :
    data_(data),
    log(log)
{
}

void FieldSamplingWriter::initialiseSampler()
{
    fieldSamplingFilter_ = vtkSmartPointer<engysResampleToImage>::New();
    scalar bounds[6];
    runTimeVis::boundsUtils::scalarArrayFromBoundBox(bounds, data_.objectData.samplingBounds);
    fieldSamplingFilter_->SetSamplingBounds(bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

    fieldSamplingFilter_->SetSamplingDimensions(
        data_.objectData.elementsX,
        data_.objectData.elementsY,
        data_.objectData.elementsZ
    );
    fieldSamplingFilter_->UseInputBoundsOff();
    fieldSamplingFilter_->SetController(vtkMultiProcessController::GetGlobalController());

    for (const word& field : data_.desiredFields)
    {
      fieldSamplingFilter_->AddField(field.c_str());
    }

    switch(data_.objectData.toleranceType.getValue())
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
    fieldSamplingFilter_->SetTolerance(data_.objectData.toleranceValue);

    pvdWriter_ = vtkSmartPointer<engysPVDWriter>::New();
    pvdWriter_->SetFileName(data_.pvdFileName.c_str());
    pvdWriter_->SetCompressorTypeToLZ4();
    pvdWriter_->SetController(vtkMultiProcessController::GetGlobalController());
    pvdWriter_->SetInputConnection(fieldSamplingFilter_->GetOutputPort());
    pvdWriter_->EncodeAppendedDataOff();
    pvdWriter_->SetWriteMultipleFilesInParallel(true);
}

void FieldSamplingWriter::write(
    const Time &currentTime,
    const runTimeVis::PostDictObjectProviderDatabase &database
)
{
    ADD_MEASURE_POINT("Starting FS write")
    if (!fieldSamplingFilter_)
    {
        initialiseSampler();
        ADD_MEASURE_POINT("Done initializing FS")
    }

    vtkDataSet *dataset = database.getDataSetForBaseItem(data_.source, currentTime.timeIndex(), currentTime.value());
    ADD_MEASURE_POINT("Got dataset")

    fieldSamplingFilter_->SetInputDataObject(dataset);
    //fieldSamplingFilter_->Update();
    //ADD_MEASURE_POINT("Updated filter")

    if (pvdWriter_)
    {
        pvdWriter_->SetCurrentTimeName(currentTime.timeName().c_str());
        pvdWriter_->Write();
        ADD_MEASURE_POINT("Wrote pvd")

        unsigned long errorCode = pvdWriter_->GetErrorCode();
        if (errorCode == vtkErrorCode::NoError)
        {
            if (log) Info << "Wrote the data to " << data_.pvdFileName << endl;
        }
        else
        {
            FatalErrorInFunction << "Error writing data to " << data_.pvdFileName << ": "
                                 << vtkErrorCode::GetStringFromErrorCode(errorCode) << abort(FatalError);
        }
    }
}

void FieldSamplingWriter::clearSampler()
{
    fieldSamplingFilter_ = nullptr;
    pvdWriter_ = nullptr;
}

} // End namespace Foam


// ************************************************************************* //
