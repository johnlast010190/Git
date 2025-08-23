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
    (c) 2020-2023 Engys Ltd.

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "runTimeVisualisationInfo.H"

#include "Utils/ParallelUtils.H"

namespace Foam::functionObjects::runTimeVis
{

RunTimeVisualisationInfo::RunTimeVisualisationInfo
(
    const word& name,
    const Time& runTime,
    const dictionary& foDict,
    const dictionary& rtppDict
)
:
    functionName_(name),
    outputFolder_( (runTimeVis::ParallelUtils::isRunningInParallel() ? runTime.path()/".." : runTime.path()) /"postProcessing"/name),
    colourMaps_(),
    renderInfo_(foDict),
    pvdInfo_(foDict),
    defaultColourLookupTablesInfo_(rtppDict),
    activeScenes_(foDict.lookup<List<word>>("activeScenes"))
{}

} // End namespace Foam