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
#include "runTimeVisualisation.H"

#include "db/Time/Time.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "signals/sigFpe.H"

#if BASIC_PROFILE == 1
#include "Utils/basicProfiler.H"
static int rtppExecutions = 0;
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

// VTK includes
#include "vtkDummyController.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam::functionObjects
{
    defineTypeNameAndDebug(runTimeVisualisation, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        runTimeVisualisation,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::runTimeVisualisation::beforeRendering()
{
    vtkMultiProcessController* controller = database_.getController();
    if (log) Info <<  "Visualising " << name() << " on "
                  << (controller ? controller->GetNumberOfProcesses() : 1) << " process"
                  << (controller ? "es" : "") << endl;

    turnFloatingPointExceptionsOff();

    database_.updateMeshDomainForTimestep(time().timeIndex(), time().value());
    database_.createItemProviders(time());
    database_.updateExternalDomainForTimestep(time().timeIndex());

#if KEEP_SCENES_IN_MEMORY
    createScenes();
#endif
}

#if KEEP_SCENES_IN_MEMORY
void Foam::functionObjects::runTimeVisualisation::createScenes() {
    if (!scenes_.empty()) {
        return;
    }

    for (const std::string& activeScene : rtppInfo_.getActiveScenes())
    {
        scenes_.emplace_back(activeScene,
                             database_,
                             rtppInfo_);
    }
}

void Foam::functionObjects::runTimeVisualisation::renderScenes()
{
    for (runTimeVis::Scene &scene : scenes_)
    {
#if BASIC_PROFILE == 1
        runTimeVis::basicProfiler::startGroup("rendering scene " + scene.getName());
#endif
        Info<< indent << "rendering scene " << scene.getName() << incrIndent << endl;
        scene.renderAndExport(time(), time_);
#if BASIC_PROFILE == 1
        runTimeVis::basicProfiler::addMeasurePoint("Finished rendering scene " + scene.getName());
        runTimeVis::basicProfiler::endGroup();
#endif
    }
}

#else
void Foam::functionObjects::runTimeVisualisation::renderScenes()
{
    for (const std::string& activeScene : rtppInfo_.getActiveScenes())
    {
        runTimeVis::Scene scene (activeScene, database_, rtppInfo_);
#if BASIC_PROFILE == 1
        runTimeVis::basicProfiler::startGroup("rendering scene " + scene.getName());
#endif
        Info<< indent << "rendering scene " << scene.getName() << incrIndent << endl;
        scene.renderAndExport(time(), time_);
#if BASIC_PROFILE == 1
        runTimeVis::basicProfiler::addMeasurePoint("Finished rendering scene " + scene.getName());
        runTimeVis::basicProfiler::endGroup();
#endif
    }
}
#endif

void Foam::functionObjects::runTimeVisualisation::afterRendering()
{
#if !KEEP_PROVIDERS_IN_MEMORY
    database_.deleteItemProviders(time());
#endif

    turnFloatingPointExceptionsOn();
}

void Foam::functionObjects::runTimeVisualisation::turnFloatingPointExceptionsOff()
{
    // Disable any floating point trapping
    // (OSMesa VTK causes floating point errors)
    sigFpe::unset();
}

void Foam::functionObjects::runTimeVisualisation::turnFloatingPointExceptionsOn()
{
    sigFpe::set();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeVisualisation::runTimeVisualisation
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    rtppFunctionObject(name, runTime, dict),
    dictionaries_(runTime, dict),
    rtppInfo_(name, runTime, dict, dictionaries_.getRtppDict()),
    database_(name, runTime, dictionaries_, debug)
{
    Info<< "Initialising runtime visualisation " << name << "..." << endl;

#if BASIC_PROFILE == 1
    runTimeVis::basicProfiler::startGroup("to initialise RTPP");
#endif

    for (const std::string& activeScene : rtppInfo_.getActiveScenes())
    {
        database_.addToRequiredScenes(activeScene);
    }

    database_.addToItemRequirements(rtppInfo_.getPvdInfo().getItemRequirements());

    Info<< "Runtime visualisation initialised" << endl;

    finishInitialization();
#if BASIC_PROFILE == 1
    runTimeVis::basicProfiler::addMeasurePoint("Finished initialising");
    runTimeVis::basicProfiler::endGroup();
    runTimeVis::basicProfiler::startGroup("initialisation of core and first iteration");
#endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::runTimeVisualisation::write()
{
#if BASIC_PROFILE == 1
    runTimeVis::basicProfiler::addMeasurePoint("Finished iterations");
    runTimeVis::basicProfiler::endGroup();
    runTimeVis::basicProfiler::startGroup("write number " + std::to_string(rtppExecutions));
    rtppExecutions++;
    runTimeVis::basicProfiler::startGroup("pre rendering operations");
#endif
    beforeRendering();
#if BASIC_PROFILE == 1
    runTimeVis::basicProfiler::addMeasurePoint("Finished pre rendering operations");
    runTimeVis::basicProfiler::endGroup();
    runTimeVis::basicProfiler::startGroup("rendering operations");
#endif
    renderScenes();
#if BASIC_PROFILE == 1
    runTimeVis::basicProfiler::addMeasurePoint("Finished rendering operations");
    runTimeVis::basicProfiler::endGroup();
    runTimeVis::basicProfiler::startGroup("post rendering operations");
#endif
    afterRendering();
#if BASIC_PROFILE == 1
    runTimeVis::basicProfiler::addMeasurePoint("Finished post rendering operations");
    runTimeVis::basicProfiler::endGroup();
    runTimeVis::basicProfiler::endGroup();
    runTimeVis::basicProfiler::addMeasurePoint("Started profile print");
    runTimeVis::basicProfiler::printResults();
    runTimeVis::basicProfiler::startGroup("core iterations");
#endif
    return true;
}

// ************************************************************************* //
