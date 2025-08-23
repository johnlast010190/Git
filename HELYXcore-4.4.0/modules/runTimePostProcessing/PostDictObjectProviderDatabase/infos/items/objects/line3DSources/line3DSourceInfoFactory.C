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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "line3DSourceInfoFactory.H"
#include "postDict/postDictKeys.H"
#include "types/line3DSourceType.H"

#include "featureLine3DSourceInfo.H"
#include "fromFileLine3DSourceInfo.H"
#include "surfaceIntersectionLine3DSourceInfo.H"
#include "straightLine3DSourceInfo.H"

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * *  Public Member Functions  * * * * * * * * * * * //

Line3DSourceInfo* Line3DSourceInfoFactory::getFromDict
(
    const dictionary& objectDict,
    const string& caseFolder
)
{
    const dictionary& sourceDict = objectDict.subDict(line3DObjectKeys::SOURCE_TYPE_DATA_KEY);
    const auto sourceType = sourceDict.lookup<Line3DSourceType>(line3DObjectKeys::TYPE_KEY);
    switch (sourceType.getValue())
    {
        case Line3DSourceType::STRAIGHT_LINE:
            return new StraightLine3DSourceInfo(sourceDict);
        case Line3DSourceType::SURFACE_INTERSECTION:
            return new SurfaceIntersectionLine3DSourceInfo(sourceDict, caseFolder);
        case Line3DSourceType::FEATURE_LINE:
            return new FeatureLine3DSourceInfo(sourceDict);
        case Line3DSourceType::FROM_FILE:
            return new FromFileLine3DSourceInfo(sourceDict, caseFolder);
        default:
            FatalError << "unknown line 3D source type " << sourceType << exit(FatalError);
            return nullptr;
    }
}

// ************************************************************************* //
}
