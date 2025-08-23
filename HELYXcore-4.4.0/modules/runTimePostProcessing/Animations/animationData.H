/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : Dev
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
    Foam::functionObjects::surfaceStatistics::AnimationData

Description
    Contains the common Animation data read from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "controlDict/controlDictKeysAN.H"
#include "types/animationType.H"
#include "types/timingMode.H"

#include "PostDictObjectProviderDatabase/baseClasses/id.H"
#include "db/dictionary/dictionary.H"
#include "baseClasses/foamField.H"
#include "meshes/primitiveShapes/point/point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
namespace Foam::functionObjects::animations
{

struct AnimationData
{
    double frameRate;
    double frameDuration;
    int width;
    int height;
    List<AnimationType> outputFormats;
    TimingMode timingMode;
    double timeStart;
    double timeEnd;
    double timeScale;
    List<word> RTPPInputs;
    string ffmpegPath;

    explicit AnimationData(const dictionary &animationDict)
    {
        frameRate = animationDict.lookupOrDefault<double>(controlKeys::FRAME_RATE_KEY, 30);
        width = animationDict.lookupOrDefault<int>(controlKeys::WIDTH_KEY, -1);
        height = animationDict.lookupOrDefault<int>(controlKeys::HEIGHT_KEY, -1);
        outputFormats = animationDict.lookup<List<AnimationType>>(controlKeys::OUTPUT_FORMAT_KEY);
        timingMode = animationDict.lookup<TimingMode>(controlKeys::TIMING_MODE_KEY);
        timeStart = animationDict.lookupOrDefault<double>(controlKeys::TIME_START_KEY, -1);
        timeEnd = animationDict.lookupOrDefault<double>(controlKeys::TIME_END_KEY, -1);
        RTPPInputs = animationDict.lookup<List<word>>(controlKeys::RTPP_INPUTS_KEY);
        ffmpegPath = animationDict.lookupOrDefault<string>(controlKeys::FFMPEG_PATH_KEY, "");
        timeScale = animationDict.lookupOrDefault<double>(controlKeys::TIME_SCALE_KEY, 0.5);
        frameDuration = animationDict.lookupOrDefault<double>(controlKeys::FRAME_DURATION_KEY, 0.5);
    }

    bool operator==(const AnimationData &other) const
    {
        return frameRate == other.frameRate &&
            width == other.width &&
            height == other.height &&
            timingMode == other.timingMode &&
            timeStart == other.timeStart &&
            timeEnd == other.timeEnd &&
            timeScale == other.timeScale &&
            frameDuration == other.frameDuration &&
            outputFormats == other.outputFormats &&
            RTPPInputs == other.RTPPInputs;
    }
};

} // End namespace Foam

// ************************************************************************* //
