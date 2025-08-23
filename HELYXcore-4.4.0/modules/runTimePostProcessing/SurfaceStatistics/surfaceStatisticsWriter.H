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

Class
    Foam::functionObjects::surfaceStatistics::SurfaceStatisticsWriter

Description
    Calculates and writes the statistical values for a given dataset

SourceFiles
    surfaceStatisticsWriter.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "primitives/Scalar/scalar/scalar.H"
#include "primitives/strings/fileName/fileName.H"
#include "baseClasses/foamField.H"
#include "vtkSmartPointer.h"
#include "db/Time/Time.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "db/functionObjects/stateFunctionObject/stateFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class engysSurfaceIntegral;
class engysCSVExporter;
class vtkAppendPolyData;
class vtkMultiProcessController;

namespace Foam::functionObjects
{

namespace runTimeVis
{
class PostDictObjectProviderDatabase;
class ReferenceFrames;
}

namespace surfaceStat
{

struct SurfaceIntegralData;

/*---------------------------------------------------------------------------*\
                  Class SurfaceStatisticsWriter Declaration
\*---------------------------------------------------------------------------*/

class SurfaceStatisticsWriter
{
private:

    vtkSmartPointer<engysSurfaceIntegral> integrator_;
    vtkSmartPointer<engysCSVExporter> csvExporter_;
    vtkSmartPointer<vtkAppendPolyData> appender_;

    const word name_;
    const SurfaceIntegralData &data_;

    class InternalWriter
    {
    public:
        InternalWriter(
            writeFile &fileWriter_,
            stateFunctionObject &fnObjectState_,
            vtkMultiProcessController *controller_
        );

        void initialize(const List<runTimeVis::foamField>& fields) const;

        void writeTime();
        void writeValueAndPosition(const word& entryName, scalar value, const vector& position);
        void writeValue(const word& entryName, scalar value);
        void writeInvalidValue(const word &entryName);

        void finalize() const;

        void updateController(vtkMultiProcessController* controller);
        [[nodiscard]] const char* getFileName() const;
        [[nodiscard]] fileName getTimeDir() const;
        [[nodiscard]] label getPrecision() const;

    private:
        writeFile &fileWriter_;
        stateFunctionObject &fnObjectState_;
        vtkMultiProcessController* controller_;
    } internalWriter_;

    const Switch& log;

    void addSourceDataToAppender(
        const runTimeVis::PostDictObjectProviderDatabase &database,
        label timeIndex,
        scalar currentTime
    );

    void writeFiles();
    void writeField(const string &fieldName);

    void initialiseIntegrator();
    void initialiseCsvExporter(const runTimeVis::ReferenceFrames &referenceFrames);
public:
    explicit SurfaceStatisticsWriter(
        const word &name,
        const SurfaceIntegralData &data,
        writeFile &fileWriter,
        stateFunctionObject& fnObjectState,
        const Switch &log
    );

    void write(
        const Time &currentTime,
        const runTimeVis::PostDictObjectProviderDatabase &database
    );

    void clearIntegrator();

    ~SurfaceStatisticsWriter() = default;
    static std::string vecToString(const vector& vector);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceStat
} // End namespace Foam::functionObjects

// ************************************************************************* //
