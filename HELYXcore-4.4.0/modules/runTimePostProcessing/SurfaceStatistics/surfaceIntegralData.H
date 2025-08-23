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
    Foam::functionObjects::surfaceStatistics::SurfaceIntegralData

Description
    Contains the common Surface Statistics data read from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "types/weighting.H"
#include "controlDict/controlDictKeys.H"

#include "PostDictObjectProviderDatabase/baseClasses/id.H"
#include "db/dictionary/dictionary.H"
#include "baseClasses/foamField.H"
#include "meshes/primitiveShapes/point/point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
namespace Foam::functionObjects::surfaceStat
{

struct SurfaceIntegralData
{
    bool flipNormal;
    Weighting weighting;
    List<runTimeVis::foamField> fields;
    List<runTimeVis::Id> sources;
    bool orientToPoint;
    point orientationPoint;
    word referenceFrame;
    bool exportResults;


    explicit SurfaceIntegralData(const dictionary &surfaceStatDict)
    {
        flipNormal = surfaceStatDict.lookupOrDefault(controlKeys::FLIP_NORMAL_KEY, false);
        weighting = surfaceStatDict.lookupOrDefault(controlKeys::WEIGHTING_KEY, Weighting());

        auto fieldWords = surfaceStatDict.lookup<List<word>>(controlKeys::FIELDS_KEY);
        for (const word &fieldWord: fieldWords)
        {
            fields.append(runTimeVis::foamField(fieldWord));
        }

        auto sourcesDictList = surfaceStatDict.lookup<List<dictionary>>(controlKeys::SOURCES_KEY);
        sources.clear();
        for (const dictionary &sourceDict: sourcesDictList)
        {
            sources.append(runTimeVis::Id(sourceDict));
        }

        orientToPoint = surfaceStatDict.found(controlKeys::ORIENTATION_POINT_KEY);
        if (orientToPoint)
        {
            orientationPoint = surfaceStatDict.lookup<point>(controlKeys::ORIENTATION_POINT_KEY);
        }

        exportResults = surfaceStatDict.lookupOrDefault<bool>(controlKeys::EXPORT_RESULTS_KEY, false);
        if (exportResults)
        {
            referenceFrame = surfaceStatDict.lookupOrDefault<word>(controlKeys::REFERENCE_FRAME_KEY, "");
        }
    }

    bool operator==(const SurfaceIntegralData &other) const
    {
        if (fields.size() != other.fields.size())
        {
            return false;
        }
        bool res = flipNormal == other.flipNormal &&
                   weighting == other.weighting &&
                   exportResults == other.exportResults;
        if (!res)
        {
            return false;
        }
        if (exportResults)
        {
            if (referenceFrame != other.referenceFrame)
            {
                return false;
            }
        }
        return std::all_of(
            fields.begin(), fields.end(),
            [&other](const runTimeVis::foamField &v) {
                return other.fields.found(v);
            }
        );
    }
};

} // End namespace Foam

// ************************************************************************* //
