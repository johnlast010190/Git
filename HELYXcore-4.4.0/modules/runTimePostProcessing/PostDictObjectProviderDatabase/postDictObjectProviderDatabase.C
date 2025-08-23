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

#include <memory>

#include "vtkMultiProcessController.h"

#include "postDictObjectProviderDatabase.H"
#include "databaseInstance/postDictObjectProviderDatabaseInstance.H"
#include "dictionaries/dictionaries.H"

#include "include/OSspecific.H"

#include "Utils/memoryMonitor.H"

#include "Utils/ParallelUtils.H"
#include "vtkMPIController.h"

#undef Log
#include "vtkLogger.h"

namespace Foam::functionObjects::runTimeVis
{

std::unique_ptr<PostDictObjectProviderDatabaseInstance> PostDictObjectProviderDatabase::instance_ = nullptr;
unsigned int PostDictObjectProviderDatabase::instanceReferenceCount_ = 0;
bool PostDictObjectProviderDatabase::createdProviders_ = false;
std::unordered_set<std::string> PostDictObjectProviderDatabase::registeredFunctionObjects_;
PostDictObjectProviderDatabase::LoggerState PostDictObjectProviderDatabase::loggerInitialized_ = NOT_INITIALIZED;
bool PostDictObjectProviderDatabase::debug_ = false;
unsigned int PostDictObjectProviderDatabase::staticReferenceCount_ = 0;

PostDictObjectProviderDatabase::PostDictObjectProviderDatabase
(
    const word &name,
    const Time &runTime,
    const Dictionaries &dictionaries,
    bool debug
)
{
    this->name = name;
    registeredFunctionObjects_.insert(name);

    if (instance_ == nullptr)
    {
        instance_ = std::make_unique<PostDictObjectProviderDatabaseInstance>(runTime, dictionaries);
        createdProviders_ = false;
    }

    localInstance_ = instance_.get();
    instanceReferenceCount_++;
}

void PostDictObjectProviderDatabase::initialiseMpi()
{
    if (ParallelUtils::isRunningInParallel())
    {
        if (!vtkMultiProcessController::GetGlobalController())
        {
            vtkMPIController* controller = vtkMPIController::New();
            // This is initialised externally by Pstream
            if (!controller->GetCommunicator())
            {
                controller->Initialize(nullptr, nullptr, 1);
            }
            vtkMultiProcessController::SetGlobalController(controller);
        }
    }
    else
    {
        vtkMultiProcessController::SetGlobalController(nullptr);
    }
}

void PostDictObjectProviderDatabase::finalizeStaticCommon()
{
    staticReferenceCount_--;
    if (staticReferenceCount_ == 0)
    {
        while (vtkMultiProcessController *controller = vtkMultiProcessController::GetGlobalController())
        {
            controller->Finalize(1);
            controller->Delete();
            vtkMultiProcessController::SetGlobalController(nullptr);
        }
    }
}

void PostDictObjectProviderDatabase::initializeStaticCommon(bool debug)
{
    staticReferenceCount_++;
    debug_ |= debug;
    initialiseMpi();
    initializeLogger(debug);
}

void PostDictObjectProviderDatabase::initializeLogger(bool debug)
{
    if (loggerInitialized_ == NOT_INITIALIZED)
    {
        vtkLogger::Init();
        vtkLogger::SetStderrVerbosity(vtkLogger::VERBOSITY_WARNING);

        if (!debug)
        {
            vtkMultiProcessController *controller = vtkMultiProcessController::GetGlobalController();
            if (!controller || controller->GetLocalProcessId() == 0)
            {
                // Only log INFO, WARNING, ERROR to "evtk.log":
                std::string fileName;
                fileName = "log/evtk.log";
                mkDir("log");
                vtkLogger::LogToFile(fileName.c_str(), vtkLogger::TRUNCATE, vtkLogger::VERBOSITY_WARNING);
            }
        }
        loggerInitialized_ = INITIALIZED;
    }

    if (debug && loggerInitialized_ == INITIALIZED)
    {
        vtkMultiProcessController *controller = vtkMultiProcessController::GetGlobalController();
        std::string fileName;
        if (controller && controller->GetNumberOfProcesses() > 1)
        {
            fileName = "log/vtkLog_" + std::to_string(controller->GetLocalProcessId()) + ".log";
        }
        else
        {
            fileName = "log/vtkLog.log";
        }
        mkDir("log");
        vtkLogger::LogToFile(fileName.c_str(), vtkLogger::TRUNCATE, vtkLogger::VERBOSITY_INFO);
        loggerInitialized_ = DEBUG;
    }
}

PostDictObjectProviderDatabase::~PostDictObjectProviderDatabase()
{
    instanceReferenceCount_--;
    if (instanceReferenceCount_ == 0)
    {
        instance_.reset(nullptr);
        createdProviders_ = false;
    }
}

void PostDictObjectProviderDatabase::addToRequiredItems(const Id &id)
{
    localInstance_->addToRequiredItems(id);
}

void PostDictObjectProviderDatabase::addToRequiredScenes(const std::string& sceneName)
{
    localInstance_->addToRequiredScenes(sceneName);
}

void PostDictObjectProviderDatabase::addToItemRequirements(const ItemRequirements &itemRequirements)
{
    localInstance_->addToItemRequirements(itemRequirements);
}

size_t PostDictObjectProviderDatabase::getSceneCount()
{
    return localInstance_->getSceneCount();
}

const SceneInfo &PostDictObjectProviderDatabase::getSceneInfo(int index) const
{
    return localInstance_->getSceneInfo(index);
}

const SceneInfo &PostDictObjectProviderDatabase::getSceneInfo(const std::string &sceneName) const
{
    return localInstance_->getSceneInfo(sceneName);
}

void PostDictObjectProviderDatabase::updateMeshDomainForTimestep(label timeIndex, scalar currentTime)
{
    localInstance_->updateMeshDomainForTimestep(timeIndex, currentTime);
}

void PostDictObjectProviderDatabase::updateExternalDomainForTimestep(label timeIndex)
{
    localInstance_->updateExternalDomainForTimestep(timeIndex);
}

vtkSmartPointer<vtkDataSet> PostDictObjectProviderDatabase::getDataSetForSceneItem(
    const ItemInfo *info,
    const std::string &sceneName,
    label timeIndex,
    scalar currentTime
) const
{
    ItemDataSetProvider* provider = localInstance_->getDataSetProviderForSceneItem(info, sceneName);
    provider->updateIfNecessary(timeIndex, currentTime);
    return provider->getDataSetOutput();
}

vtkSmartPointer<vtkPolyData> PostDictObjectProviderDatabase::getOutlineForSceneItem(
    const ItemInfo *info,
    const std::string &sceneName,
    label timeIndex,
    scalar currentTime
) const
{
    ItemDataSetProvider* provider = localInstance_->getDataSetProviderForSceneItem(info, sceneName);
    provider->updateIfNecessary(timeIndex, currentTime);
    return provider->getDataSetOutline();
}

vtkSmartPointer<vtkDataSet> PostDictObjectProviderDatabase::getDataSetForBaseItem(
    const Id &id,
    label timeIndex,
    scalar currentTime
) const
{
    ItemDataSetProvider* provider = localInstance_->getDataSetProviderForBaseItem(id);
    provider->updateIfNecessary(timeIndex, currentTime);
    return provider->getDataSetOutput();
}

const ColourLookupTablesInfo& PostDictObjectProviderDatabase::getBaseColorLookupTablesInfo() const
{
    return localInstance_->getBaseColorLookupTablesInfo();
}

const FoamMeshes& PostDictObjectProviderDatabase::getFoamMeshes() const
{
    return localInstance_->getFoamMeshes();
}

const ExternalFields& PostDictObjectProviderDatabase::getExternalFields() const
{
    return localInstance_->getExternalFields();
}

vtkMultiProcessController* PostDictObjectProviderDatabase::getController() const
{
    return localInstance_->getController();
}

void PostDictObjectProviderDatabase::createItemProviders(const Time& vtkNotUsed(time))
{
    if (!createdProviders_)
    {
        if (debug_) vtkLogF(INFO, "Creating item providers. Memory consumption: current: %zu mb; peak: %zu mb", MemoryMonitor::getCurrentRSSMByte(), MemoryMonitor::getPeakRSSMByte());
        localInstance_->createItemProviders();
        createdProviders_ = true;
        if (debug_) vtkLogF(INFO, "Created item providers. Memory consumption: current: %zu mb; peak: %zu mb", MemoryMonitor::getCurrentRSSMByte(), MemoryMonitor::getPeakRSSMByte());
    }
}

void PostDictObjectProviderDatabase::deleteItemProviders(const Time& time)
{
    if (isLastFunctionObject(name, time))
    {
        localInstance_->deleteItemProviders();
        createdProviders_ = false;
        if (debug_) vtkLogF(INFO, "Deleted item providers. Memory consumption: current: %zu mb; peak: %zu mb", MemoryMonitor::getCurrentRSSMByte(), MemoryMonitor::getPeakRSSMByte());
    }
}

bool PostDictObjectProviderDatabase::isLastFunctionObject(const word& name, const Time &time)
{
    const functionObjectList& foList = time.functionObjects();

    word lastFO;
    forAll(foList, n)
    {
        const word& foName = foList[n].name();
        if (registeredFunctionObjects_.find(foName) != registeredFunctionObjects_.end())
        {
            lastFO = foName;
        }
    }
    return lastFO == name;
}

const ReferenceFrames& PostDictObjectProviderDatabase::getReferenceFrames() const
{
    return localInstance_->getReferenceFrames();
}

}

// ************************************************************************* //
