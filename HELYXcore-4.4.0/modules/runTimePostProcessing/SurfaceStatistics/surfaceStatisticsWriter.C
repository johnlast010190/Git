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
    (c) 2022-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "surfaceStatisticsWriter.H"

#include "surfaceIntegralData.H"
#include "compileOptions.H"

#include "Utils/vtkDataHandlingTools.H"
#include "storage/referenceFrames.H"
#include "postDictObjectProviderDatabase.H"

#include "engysSurfaceIntegral.h"
#include "engysCSVExporter.h"
#include "vtkAppendPolyData.h"
#include "vtkMultiProcessController.h"

#include "include/OSspecific.H"

namespace Foam::functionObjects::surfaceStat
{

SurfaceStatisticsWriter::InternalWriter::InternalWriter(
    writeFile &fileWriter_,
    stateFunctionObject &fnObjectState_,
    vtkMultiProcessController *controller_
) : fileWriter_(fileWriter_), fnObjectState_(fnObjectState_), controller_(controller_)
{}

void SurfaceStatisticsWriter::InternalWriter::initialize(const List<runTimeVis::foamField>& fields) const
{
    // .dat file header
    if (!controller_ || controller_->GetLocalProcessId() == 0)
    {
        fileWriter_.writeCommented(fileWriter_.file(), "Time");

        fileWriter_.writeDelimited(fileWriter_.file(), "flux");

        for (const string& field: fields)
        {
            // TODO if it's a vector, there should be components
            fileWriter_.writeDelimited(fileWriter_.file(), field + "min");
            fileWriter_.writeDelimited(fileWriter_.file(), "at location");
            fileWriter_.writeDelimited(fileWriter_.file(), field + "max");
            fileWriter_.writeDelimited(fileWriter_.file(), "at location");
            fileWriter_.writeDelimited(fileWriter_.file(), field + "mean");
            fileWriter_.writeDelimited(fileWriter_.file(), field + "stDev");
        }

        fileWriter_.file() << endl;
    }
}

void SurfaceStatisticsWriter::InternalWriter::writeTime()
{
    if (!controller_ || controller_->GetLocalProcessId() == 0)
    {
        fileWriter_.writeTime(fileWriter_.file());
    }
}

void SurfaceStatisticsWriter::InternalWriter::writeValueAndPosition(const word& entryName, scalar value, const vector& position)
{

    writeValue(entryName, value);

    if (!controller_ || controller_->GetLocalProcessId() == 0)
    {
        Ostream &os = fileWriter_.file();
        fileWriter_.writeDelimited(os, vecToString(position));
    }

    fnObjectState_.setResult
        (
            entryName+"_position",
            position
        );
}

void SurfaceStatisticsWriter::InternalWriter::writeValue(const word& entryName, scalar value)
{
    if (!controller_ || controller_->GetLocalProcessId() == 0)
    {
        Ostream &os = fileWriter_.file();
        fileWriter_.writeDelimited(os, value);
    }

    fnObjectState_.setResult(entryName, value);
}

void SurfaceStatisticsWriter::InternalWriter::writeInvalidValue(const word& entryName)
{
    Ostream &os = fileWriter_.file();
    fileWriter_.writeDelimited(os, "NegativeWeights");

    fnObjectState_.setResult(entryName, scalar(-1));
}

void SurfaceStatisticsWriter::InternalWriter::finalize() const
{
    if (!controller_ || controller_->GetLocalProcessId() == 0)
    {
        fileWriter_.file() << endl;
    }
}

void SurfaceStatisticsWriter::InternalWriter::updateController(vtkMultiProcessController* controller)
{
    this->controller_ = controller;
}

const char* SurfaceStatisticsWriter::InternalWriter::getFileName() const
{
    return fileWriter_.file().name().c_str();
}

label SurfaceStatisticsWriter::InternalWriter::getPrecision() const
{
    return fileWriter_.charWidth() - Foam::functionObjects::writeFile::addChars;
}

fileName SurfaceStatisticsWriter::InternalWriter::getTimeDir() const
{
    return fileWriter_.baseTimeDir();
}

// ---------------------------------------------------------------------
SurfaceStatisticsWriter::SurfaceStatisticsWriter(
    const word &name,
    const SurfaceIntegralData &data,
    writeFile &fileWriter,
    stateFunctionObject& fnObjectState,
    const Switch &log
)
    :
    name_(name),
    data_(data),
    internalWriter_(fileWriter, fnObjectState, vtkMultiProcessController::GetGlobalController()),
    log(log)
{
    if (log) Info << "    Logging surface statistics to file: " << internalWriter_.getFileName() << endl;
    internalWriter_.initialize(data_.fields);
}

void SurfaceStatisticsWriter::initialiseIntegrator()
{
    integrator_ = vtkSmartPointer<engysSurfaceIntegral>::New();
    integrator_->SetUseFluxAsWeight(data_.weighting.IsFlux());
    integrator_->SetFlipNormal(data_.flipNormal);
    for (const string &desiredFieldString: data_.fields)
    {
        runTimeVis::foamField desiredField(desiredFieldString);
        integrator_->AddField(desiredField.c_str());
    }
    integrator_->SetOrientToPoint(data_.orientToPoint);
    if (data_.orientToPoint)
    {
        integrator_->SetOrientationPoint(
            data_.orientationPoint.x(),
            data_.orientationPoint.y(),
            data_.orientationPoint.z());
    }
}

void SurfaceStatisticsWriter::initialiseCsvExporter(const runTimeVis::ReferenceFrames& referenceFrames)
{
    if (data_.exportResults)
    {
        csvExporter_ = vtkSmartPointer<engysCSVExporter>::New();

        csvExporter_->SetController(vtkMultiProcessController::GetGlobalController());

        csvExporter_->SetPrecision(internalWriter_.getPrecision());

        csvExporter_->SetFlipNormal(data_.flipNormal);
        for (const string &desiredFieldString: data_.fields)
        {
            runTimeVis::foamField desiredField(desiredFieldString);
            csvExporter_->AddField(desiredField.c_str());
        }
        csvExporter_->SetOrientToPoint(data_.orientToPoint);
        if (data_.orientToPoint)
        {
            csvExporter_->SetOrientationPoint(
                data_.orientationPoint.x(),
                data_.orientationPoint.y(),
                data_.orientationPoint.z());
        }

        if (data_.referenceFrame.empty())
        {
            csvExporter_->UseReferenceFrameOff();
        }
        else
        {
            csvExporter_->UseReferenceFrameOn();
            const runTimeVis::ReferenceFrame &referenceFrame = referenceFrames.getReferenceFrame(data_.referenceFrame);
            if (referenceFrame.getType().getValue() != runTimeVis::CoordinateSystemType::CARTESIAN)
            {
                FatalErrorInFunction << "Invalid coordinate system type for " << data_.referenceFrame
                                     << ": only cartesian is supported." << abort(FatalError);
                return;
            }
            point origin = referenceFrame.getOrigin();
            csvExporter_->SetReferenceFrameOrigin(origin[0], origin[1], origin[2]);
            vector e1 = referenceFrame.e1();
            vector e2 = referenceFrame.e2();
            csvExporter_->SetReferenceFrameXDirection(e1[0], e1[1], e1[2]);
            csvExporter_->SetReferenceFrameYDirection(e2[0], e2[1], e2[2]);
        }
    }
}

void SurfaceStatisticsWriter::writeFiles()
{
    if (log) Info << "    Total Area: " << integrator_->GetTotalArea()
                  << "    Total Flux: " << integrator_->GetTotalCellFlux()
                  << endl;

    internalWriter_.writeTime();
    internalWriter_.writeValue("flux", static_cast<scalar>(integrator_->GetTotalCellFlux()));

    for (const string &fieldName: data_.fields)
    {
        writeField(fieldName);
    }

    internalWriter_.finalize();
}

void SurfaceStatisticsWriter::writeField(const string &fieldName)
{
    enum StandardDeviationErrorCodes
    {
        VALID,
        NEGATIVE_WEIGHTS
    };

    const static bool isPointValue = false;

    std::string actualFieldName;
    if (integrator_->IsFieldVector(fieldName.c_str(), isPointValue))
    {
        actualFieldName = fieldName + "-Mag";
    }
    else
    {
        actualFieldName = fieldName;
    }
    const char* fieldNameC = actualFieldName.c_str();

    word minName("min("+actualFieldName+")");
    double* dLocation = integrator_->GetFieldMinLocation(fieldNameC, isPointValue);
    internalWriter_.writeValueAndPosition(minName,
                          static_cast<scalar>(integrator_->GetFieldMin(fieldNameC, isPointValue)),
                          vector(
                              static_cast<scalar>(dLocation[0]),
                              static_cast<scalar>(dLocation[1]),
                              static_cast<scalar>(dLocation[2])
                              )
                          );

    word maxName("max("+actualFieldName+")");
    dLocation = integrator_->GetFieldMaxLocation(fieldNameC, isPointValue);
    internalWriter_.writeValueAndPosition(maxName,
                          static_cast<scalar>(integrator_->GetFieldMax(fieldNameC, isPointValue)),
                          vector(
                              static_cast<scalar>(dLocation[0]),
                              static_cast<scalar>(dLocation[1]),
                              static_cast<scalar>(dLocation[2])
                              )
                          );

    word meanName("mean("+actualFieldName+")");
    internalWriter_.writeValue(meanName, static_cast<scalar>(integrator_->GetFieldAverage(fieldNameC, isPointValue)));

    word stdDevName("stdDev("+actualFieldName+")");
    if (integrator_->GetStandardDeviationErrorCode() == VALID)
    {
        internalWriter_.writeValue(stdDevName, static_cast<scalar>(integrator_->GetFieldStd(fieldNameC, isPointValue)));
    }
    else
    {
        internalWriter_.writeInvalidValue(stdDevName);
    }
}

std::string SurfaceStatisticsWriter::vecToString(const vector& vector)
{
    return "(" + std::to_string(vector[0]) + " " + std::to_string(vector[1]) + " " + std::to_string(vector[2]) + ")";
}

void SurfaceStatisticsWriter::write(
    const Time &currentTime,
    const runTimeVis::PostDictObjectProviderDatabase &database
)
{
    internalWriter_.updateController(database.getController());
    if (!integrator_)
    {
        initialiseIntegrator();
    }

    if (!appender_)
    {
        appender_ = vtkSmartPointer<vtkAppendPolyData>::New();
#if KEEP_PROVIDERS_IN_MEMORY
        addSourceDataToAppender(data, database, currentTime.timeIndex());
    }
#else
    }
    addSourceDataToAppender(database, currentTime.timeIndex(), currentTime.value());
#endif

    appender_->Update();
    vtkSmartPointer<vtkPolyData> appended = appender_->GetOutput();
    appended->RemoveGhostCells();

    integrator_->SetInputData(appended);

    if (data_.exportResults)
    {
        if (!csvExporter_) {
            initialiseCsvExporter(database.getReferenceFrames());
        }
        csvExporter_->SetInputData(appended);
        fileName outputDir = internalWriter_.getTimeDir();
        mkDir(outputDir);
        fileName outputFileName = internalWriter_.getTimeDir() / name_ + ".csv";
        csvExporter_->SetFileName(outputFileName.c_str());
        csvExporter_->Write();
    }

    integrator_->Update();

    writeFiles();
}

void SurfaceStatisticsWriter::addSourceDataToAppender(
    const runTimeVis::PostDictObjectProviderDatabase &database,
    label timeIndex,
    scalar currentTime
)
{
    appender_->RemoveAllInputs();

    for (const runTimeVis::Id &source: data_.sources)
    {
        vtkSmartPointer<vtkDataSet> dataSet = database.getDataSetForBaseItem(source, timeIndex, currentTime);

        vtkSmartPointer<vtkPolyData> polyData = vtk::Tools::dataObjectToPolyData(dataSet);

        appender_->AddInputDataObject(polyData);
    }
}

void SurfaceStatisticsWriter::clearIntegrator()
{
    integrator_ = nullptr;
    appender_ = nullptr;
    csvExporter_ = nullptr;
}

} // End namespace Foam


// ************************************************************************* //
