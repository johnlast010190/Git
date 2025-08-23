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
#include "animation.H"

#include "compileOptions.H"
#include "types/timingMode.H"

#include "db/Time/Time.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/fileOperations/fileOperation/fileOperation.H"
#include "include/OSspecific.H"

#include "engysAnimationFilter.h"

#include "signals/sigFpe.H"

#if BASIC_PROFILE == 1
#include "Utils/basicProfiler.H"
static int anExecutions = 0;
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
defineTypeNameAndDebug(Animation, 0);

addToRunTimeSelectionTable
(
    functionObject,
    Animation,
    dictionary
);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::Animation::beforeWriting()
{
    if (log) Info <<  "Generating animation " << name() << endl;
    turnFloatingPointExceptionsOff();

    this->animationFilter = vtkSmartPointer<engysAnimationFilter>::New();

    animationFilter->SetFilePath((animationFilePath_ / name_).c_str());
    animationFilter->SetFrameRate(data_.frameRate);
    animationFilter->SetWidth(data_.width);
    animationFilter->SetHeight(data_.height);
    switch(data_.timingMode.getValue())
    {
        case animations::TimingMode::TIME:
            animationFilter->SetTimingModeToTimeScale();
            animationFilter->SetTimeScale(data_.timeScale);
            break;
        case animations::TimingMode::DURATION:
            animationFilter->SetTimingModeToFrameDuration();
            animationFilter->SetFrameDuration(data_.frameDuration);
            break;
    }
    animationFilter->SetTimeStart(data_.timeStart);
    animationFilter->SetTimeEnd(data_.timeEnd);
    animationFilter->SetFFMPEGFilePath(data_.ffmpegPath.c_str());

    for (const word& rtppFO : data_.RTPPInputs)
    {
        addRTPPImageFolder(rtppFO, 0);
    }
}

void Foam::functionObjects::Animation::turnFloatingPointExceptionsOff()
{
    // Disable any floating point trapping
    // (OSMesa VTK causes floating point errors)
    sigFpe::unset();
}

void Foam::functionObjects::Animation::addRTPPImageFolder(const word& rtppName, double timeOffset) const
{
    fileName rtppFolder = postProcessingFolder_ / rtppName;
    DebugInformation << "Adding rtpp folder " << rtppFolder << endl;

    const fileOperation& handler = fileHandler();
    fileNameList files = handler.readDir(rtppFolder, fileName::DIRECTORY);
    DebugInformation << "found " << files.size() << " files" << endl;
    std::sort(files.begin(), files.end());
    forAll(files, i)
    {
        std::string sceneName = files[i];
        fileName folder = rtppFolder / sceneName;
        DebugInformation << "Adding scene folder " << folder << endl;
        addRTPPSceneImageFolder(rtppName, sceneName, folder, timeOffset);
    }
}

void Foam::functionObjects::Animation::addRTPPSceneImageFolder(const word& rtppName, const std::string& sceneName, const fileName& sceneFolder, double timeOffset) const
{
    animationFilter->AddImageFolder(rtppName.c_str(), sceneName.c_str(), sceneFolder.c_str(), sceneFolder.name().c_str(), timeOffset);
}

void Foam::functionObjects::Animation::writeData()
{
    for (const auto& value : data_.outputFormats)
    {
        animationFilter->SetOutputFormat(engysAnimationFilter::AnimationOutputFormat(value.getValue()));
        fileName animationFilePath = animationFilePath_ / name_ + "." + value.getExtension();
//        rm(animationFilePath);
        DebugInformation << "Creating animation at " << animationFilePath << endl;

        animationFilter->CreateAnimation();
    }
}

void Foam::functionObjects::Animation::afterWriting()
{
    this->animationFilter = nullptr;
    turnFloatingPointExceptionsOn();
}

void Foam::functionObjects::Animation::turnFloatingPointExceptionsOn()
{
    sigFpe::set();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Animation::Animation
    (
        const word &name,
        const Time &runTime,
        const dictionary &dict
    )
    :
    rtppFunctionObject(name, runTime, dict),
    data_(dict)
{
    Animation::read(dict);
    Info<< "Initialising animation " << name << "..." << endl;
    postProcessingFolder_ = ( (Pstream::parRun() ? runTime.path()/".." : runTime.path()) /"postProcessing");
    name_ = name;

#if BASIC_PROFILE == 1
    runTimeVis::basicProfiler::startGroup("to initialise AN");
#endif

    animationFilePath_ = postProcessingFolder_ / name_;
    mkDir(animationFilePath_);

    Info<< "Animation initialised" << endl;

    finishInitialization();
#if BASIC_PROFILE == 1
    runTimeVis::basicProfiler::addMeasurePoint("Finished initialising");
    runTimeVis::basicProfiler::endGroup();
    runTimeVis::basicProfiler::startGroup("initialisation of core and first iteration");
#endif
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::Animation::write()
{
    if (Pstream::parRun() && !Pstream::master()) return true;
#if BASIC_PROFILE == 1
    runTimeVis::basicProfiler::addMeasurePoint("Finished iterations");
    runTimeVis::basicProfiler::endGroup();
    runTimeVis::basicProfiler::startGroup("write number " + std::to_string(anExecutions));
    anExecutions++;
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
