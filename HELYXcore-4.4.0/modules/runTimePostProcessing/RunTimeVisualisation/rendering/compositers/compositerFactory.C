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
    (c) 2020-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "compositerFactory.H"

#include "compileOptions.H"

#include "dummyCompositer.H"
#include "basicCompositer.H"
#include "edfCompositer.H"
#include "masterGatherCompositer.H"
#include "infos/renderInfo.H"
#include "Utils/ParallelUtils.H"

#if TRANSPARENT_COMPOSITER_TO_USE == ICE_T_COMPOSITER
#include "iceTCompositer.H"
#elif TRANSPARENT_COMPOSITER_TO_USE == ALL_GATHER_TRANSPARENT_COMPOSITER
#include "allGatherTransparentCompositer.H"
#elif TRANSPARENT_COMPOSITER_TO_USE == PARALLEL_DEPTH_PEEL_COMPOSITER
#include "parallelDepthPeelCompositer.H"
#endif

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

#define make(synchronizer) std::unique_ptr<Compositer>(synchronizer)

std::unique_ptr<Compositer> CompositerFactory::createCompositer
(
        vtkMultiProcessController* controller,
        bool hasTransparency,
        const RenderInfo& renderInfo
)
{
    //return make(ParallelDepthPeelCompositer(controller));
    if (renderInfo.exportEdf)
    {
        return make(new EDFCompositer(controller, renderInfo));
    }
    if (!ParallelUtils::isRunningInParallel())
    {
        return make(new DummyCompositer(renderInfo));
    }
    if (renderInfo.exportX3d)
    {
        return make(new MasterGatherCompositer(controller, renderInfo));
    }
    else if (hasTransparency)
    {
#if TRANSPARENT_COMPOSITER_TO_USE == ICE_T_COMPOSITER
        return make(IceTCompositer(controller, renderInfo));
#elif TRANSPARENT_COMPOSITER_TO_USE == ALL_GATHER_TRANSPARENT_COMPOSITER
        return make(AllGatherTransparentCompositer(controller, renderInfo));
#elif TRANSPARENT_COMPOSITER_TO_USE == PARALLEL_DEPTH_PEEL_COMPOSITER
        return make(new ParallelDepthPeelCompositer(controller, renderInfo));
#endif
    }
    else
    {
        return make(new BasicCompositer(controller, renderInfo));
    }
}

} // End namespace



// ************************************************************************* //
