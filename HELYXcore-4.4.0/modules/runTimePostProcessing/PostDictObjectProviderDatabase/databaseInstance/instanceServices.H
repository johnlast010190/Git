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
    (c) 2025 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::InstanceServices

Description
    Provides access to methods shared and accessible by all providers

SourceFiles
    instanceServices.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "primitives/ints/label/label.H"
#include "primitives/Scalar/scalar/scalar.H"
#include "vtkSmartPointer.h"
#include <vector>
#include <unordered_map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declarations
class vtkMultiProcessController;
class engysNewSetGhostCells;

namespace Foam::functionObjects::runTimeVis
{
class PostDictObjectProviderDatabaseInstance;
class ItemDataSetProvider;

class InstanceServicesRegion
{
public:
    InstanceServicesRegion(std::string region, PostDictObjectProviderDatabaseInstance& instance);

    ~InstanceServicesRegion() = default;

    void updateIfNecessary(label currentTimeIndex, scalar currentTime);

    void reset();

    engysNewSetGhostCells* getInitializedGhostCellsFilter();

private:
    std::string region_;
    PostDictObjectProviderDatabaseInstance& instance_;

    vtkSmartPointer<engysNewSetGhostCells> ghostCells_;

    std::unique_ptr<std::vector<ItemDataSetProvider*>> externalMeshBoundaries_;

    label latestProcessedTimeIndex_;
    scalar currentTime_;

    void initializeGhostCellsFilter();
};

/*---------------------------------------------------------------------------*\
          Class InstanceServices Declaration
\*---------------------------------------------------------------------------*/

class InstanceServices
{
public:
    explicit InstanceServices(PostDictObjectProviderDatabaseInstance& instance);

    ~InstanceServices() = default;

    void updateIfNecessary(label currentTimeIndex, scalar currentTime);

    void reset();

    InstanceServicesRegion* getInstanceServicesForRegion(const std::string& region);

//    void updateProvidersIfNecessary(const std::vector<ItemDataSetProvider*>& providers);

private:
    PostDictObjectProviderDatabaseInstance& instance_;

    std::unordered_map<std::string /*region*/, std::unique_ptr<InstanceServicesRegion>> regions_;

    label latestProcessedTimeIndex_;
    scalar currentTime_;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
