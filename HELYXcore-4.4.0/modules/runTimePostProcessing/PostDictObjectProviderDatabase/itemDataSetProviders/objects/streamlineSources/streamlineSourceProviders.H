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

Class
    Foam::functionObjects::runTimeVis::StreamlineSourceProvider

Description
    A base class for the streamline source providers, and definition of the different providers

SourceFiles
    lineSourceProvider.C
    patchSourceProvider.C
    pointCloudSourceProvider.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "vtkSmartPointer.h"
#include "dataStructs/objects/streamlineSources/lineSourceData.H"
#include "dataStructs/objects/streamlineSources/patchSourceData.H"
#include "dataStructs/objects/streamlineSources/pointCloudSourceData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations

namespace Foam::functionObjects::runTimeVis
{

class StreamlineSourceProvider
{
public:
        virtual vtkSmartPointer<vtkPolyData> getSeedPoints() = 0;

        [[nodiscard]] virtual bool acceptsPatchSource() const
        {
            return false;
        }

        virtual void setPatchSource(vtkDataObject* vtkNotUsed(patch))
        {
            // Nothing to do
        }

protected:
        StreamlineSourceProvider() = default;
        virtual ~StreamlineSourceProvider() = default;
};


// Definition of the source provider classes


class LineSourceProvider : public StreamlineSourceProvider
{
private:
    vtkSmartPointer<vtkPolyData> seedPoints_;
public:
    explicit LineSourceProvider(const LineSourceData& data);
    vtkSmartPointer<vtkPolyData> getSeedPoints() override;
};



class PointCloudSourceProvider : public StreamlineSourceProvider
{
private:
    vtkSmartPointer<vtkPolyData> seedPoints_;
public:
    explicit PointCloudSourceProvider(const PointCloudSourceData& data);
    vtkSmartPointer<vtkPolyData> getSeedPoints() override;
};



class PatchSourceProvider : public StreamlineSourceProvider
{
private:
    vtkSmartPointer<vtkPolyData> patchSource_;
    PatchSourceData data_;
public:
    explicit PatchSourceProvider(PatchSourceData data);
    vtkSmartPointer<vtkPolyData> getSeedPoints() override;
    [[nodiscard]] bool acceptsPatchSource() const override;
    void setPatchSource(vtkDataObject* patch) override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
