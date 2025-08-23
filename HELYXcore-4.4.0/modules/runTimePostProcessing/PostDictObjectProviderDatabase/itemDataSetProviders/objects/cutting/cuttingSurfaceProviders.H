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
    Foam::functionObjects::runTimeVis::CuttingSurfaceProviders

Description
    A base class for the cutting surface providers, and definition of the different providers

SourceFiles
    cuttingPlaneProvider.C
    cuttingSphereProvider.C
    cuttingCylinderProvider.C
    cuttingBoxProvider.C
    cuttingFileProvider.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "vtkSmartPointer.h"
#include "dataStructs/objects/cutting/cuttingPlaneData.H"
#include "dataStructs/objects/cutting/cuttingSphereData.H"
#include "dataStructs/objects/cutting/cuttingCylinderData.H"
#include "dataStructs/objects/cutting/cuttingBoxData.H"
#include "dataStructs/objects/cutting/cuttingFileData.H"

#include "vtkImplicitFunction.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class vtkImplicitFunction;
class vtkDataSet;
class vtkPolyData;

namespace Foam::functionObjects::runTimeVis
{

class CuttingSurfaceProvider
{
public:
        [[nodiscard]] vtkImplicitFunction* getImplicitFunction() const
        {
            return implicitFunction_;
        };

protected:
        vtkSmartPointer<vtkImplicitFunction> implicitFunction_;
        CuttingSurfaceProvider() = default;
        virtual ~CuttingSurfaceProvider() = default;
};

// Definition of the provider classes for the different cutting types
class CuttingPlaneProvider : public CuttingSurfaceProvider {
public:
    explicit CuttingPlaneProvider(const CuttingPlaneData &data);
    ~CuttingPlaneProvider() override = default;
};

class CuttingSphereProvider : public CuttingSurfaceProvider {
public:
    explicit CuttingSphereProvider(const CuttingSphereData &data);
    ~CuttingSphereProvider() override = default;
};

class CuttingCylinderProvider : public CuttingSurfaceProvider {
public:
    explicit CuttingCylinderProvider(const CuttingCylinderData &data);
    ~CuttingCylinderProvider() override = default;
};

class CuttingBoxProvider : public CuttingSurfaceProvider {
public:
    explicit CuttingBoxProvider(const CuttingBoxData &data);
    ~CuttingBoxProvider() override = default;
};

class CuttingFileProvider : public CuttingSurfaceProvider
{
public:
    explicit CuttingFileProvider(const CuttingFileData& data);
    ~CuttingFileProvider() override = default;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
