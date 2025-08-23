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
#include "fieldSampling.H"

#include "compileOptions.H"
#include "Utils/ParallelUtils.H"

#include "db/Time/Time.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "signals/sigFpe.H"

#if BASIC_PROFILE == 1
#include "Utils/basicProfiler.H"
static int fsExecutions = 0;
#endif

#if defined(WIN64) || defined(WIN32)
#include "vtkAutoInit.h"
#include "vtkRenderingOpenGLConfigure.h"
VTK_MODULE_INIT(vtkInteractionStyle)
VTK_MODULE_INIT(vtkRenderingFreeType)
VTK_MODULE_INIT(vtkFiltersParallelFlowPaths)
VTK_MODULE_INIT(vtkIOMPIImage)
// We only support OpenGL2 on Windows
// #ifdef VTK_OPENGL2
    VTK_MODULE_INIT(vtkRenderingContextOpenGL2)
    VTK_MODULE_INIT(vtkRenderingOpenGL2)
    VTK_MODULE_INIT(vtkRenderingVolumeOpenGL2)
// #else
//     VTK_MODULE_INIT(vtkRenderingContextOpenGL)
//     VTK_MODULE_INIT(vtkRenderingOpenGL)
//     VTK_MODULE_INIT(vtkRenderingVolumeOpenGL)
// #endif
#endif


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


namespace Foam::functionObjects
{
defineTypeNameAndDebug(fieldSampling, 0);

addToRunTimeSelectionTable
(
    functionObject,
    fieldSampling,
    dictionary
);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldSampling::beforeWriting()
{
    vtkMultiProcessController* controller = database_.getController();
    if (log) Info <<  "Calculating field sampling " << name() << " on "
                  << (controller ? controller->GetNumberOfProcesses() : 1) << " process"
                  << (controller ? "es" : "") << endl;

    turnFloatingPointExceptionsOff();

    database_.updateMeshDomainForTimestep(time().timeIndex(), time().value());
    database_.createItemProviders(time());
    database_.updateExternalDomainForTimestep(time().timeIndex());
}

void Foam::functionObjects::fieldSampling::writeData()
{
    writer_.write(time(), database_);
}

void Foam::functionObjects::fieldSampling::afterWriting()
{
#if !KEEP_PROVIDERS_IN_MEMORY
    database_.deleteItemProviders(time());
#endif

#if !KEEP_SCENES_IN_MEMORY
    writer_.clearSampler();
#endif

    turnFloatingPointExceptionsOn();
}

void Foam::functionObjects::fieldSampling::turnFloatingPointExceptionsOff()
{
    // Disable any floating point trapping
    // (OSMesa VTK causes floating point errors)
    sigFpe::unset();
}

void Foam::functionObjects::fieldSampling::turnFloatingPointExceptionsOn()
{
    sigFpe::set();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldSampling::fieldSampling
    (
        const word &name,
        const Time &runTime,
        const dictionary &dict
    )
    :
    rtppFunctionObject(name, runTime, dict),
    dictionaries_(runTime, dict),
    database_(name, runTime, dictionaries_, readDebugFlag(dict)),
    data_(dict, (runTimeVis::ParallelUtils::isRunningInParallel() ? runTime.path() / ".." : runTime.path()), name),
    writer_(data_, log)
{
    Info<< "Initialising field sampling " << name << "..." << endl;

#if BASIC_PROFILE == 1
    runTimeVis::basicProfiler::startGroup("to initialise FS");
#endif

    database_.addToRequiredItems(data_.source);

    runTimeVis::ItemRequirements itemRequirements;
    for (const word& field : data_.desiredFields)
    {
        itemRequirements.addToRequiredFields(runTimeVis::foamField(field));
    }
    database_.addToItemRequirements(itemRequirements);

    Info<< "Field sampling initialised" << endl;

    finishInitialization();
#if BASIC_PROFILE == 1
    runTimeVis::basicProfiler::addMeasurePoint("Finished initialising");
    runTimeVis::basicProfiler::endGroup();
    runTimeVis::basicProfiler::startGroup("initialisation of core and first iteration");
#endif
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldSampling::write()
{
#if BASIC_PROFILE == 1
    runTimeVis::basicProfiler::addMeasurePoint("Finished iterations");
    runTimeVis::basicProfiler::endGroup();
    runTimeVis::basicProfiler::startGroup("write number " + std::to_string(fsExecutions));
    fsExecutions++;
    runTimeVis::basicProfiler::startGroup("pre writing operations");
#endif
    beforeWriting();
#if BASIC_PROFILE == 1
    runTimeVis::basicProfiler::addMeasurePoint("Finished pre writing operations");
    runTimeVis::basicProfiler::endGroup();
    runTimeVis::basicProfiler::startGroup("writing operations");
#endif
    writeData();
#if BASIC_PROFILE == 1
    runTimeVis::basicProfiler::addMeasurePoint("Finished writing operations");
    runTimeVis::basicProfiler::endGroup();
    runTimeVis::basicProfiler::startGroup("post writing operations");
#endif
    afterWriting();
#if BASIC_PROFILE == 1
    runTimeVis::basicProfiler::addMeasurePoint("Finished post writing operations");
    runTimeVis::basicProfiler::endGroup();
    runTimeVis::basicProfiler::endGroup();
    runTimeVis::basicProfiler::addMeasurePoint("Started profile print");
    runTimeVis::basicProfiler::printResults();
    runTimeVis::basicProfiler::startGroup("core iterations");
#endif
    return true;
}

// ************************************************************************* //
