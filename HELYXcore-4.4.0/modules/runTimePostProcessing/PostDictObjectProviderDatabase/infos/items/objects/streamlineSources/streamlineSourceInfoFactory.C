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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "streamlineSourceInfoFactory.H"

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * *  Public Member Functions  * * * * * * * * * * * //

TemplatedStreamlineSourceInfo* StreamlineSourceInfoFactory::getFromDict
(
    const dictionary& streamlinesDict
)
{
    const auto type = streamlinesDict
            .subDict(streamlinesSourceKeys::SOURCE_DICT_KEY)
            .lookup<StreamlineSourceType>(streamlinesSourceKeys::TYPE_KEY);
    switch(type.getValue())
    {
    case StreamlineSourceType::LINE_SOURCE:
        return new LineSourceInfo(streamlinesDict);
    case StreamlineSourceType::PATCH_SOURCE:
        return new PatchSourceInfo(streamlinesDict);
    case StreamlineSourceType::POINT_CLOUD_SOURCE:
        return new PointCloudSourceInfo(streamlinesDict);
    }
    // Error, even though the code should never get here
    FatalError << "Unknown streamline source type: " << type << exit(FatalError);
    return nullptr;
}

// ************************************************************************* //
}
