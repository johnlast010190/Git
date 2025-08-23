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

#include "surfaceIntersectionLine3DSourceInfo.H"
#include "itemDataSetProviders/objects/line3d/surfaceIntersectionLine3DDataSetProvider.H"
#include "vtkImplicitFunction.h"
#include "postDict/postDictKeys.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

SurfaceIntersectionLine3DSourceInfo::SurfaceIntersectionLine3DSourceInfo(const dictionary& line3DDict, const string& caseFolder) :
cuttingInfo_(line3DDict.subDict(line3DObjectKeys::SURFACE_INTERSECTION_KEY), caseFolder)
{}

void SurfaceIntersectionLine3DSourceInfo::readDict(const dictionary& line3DDict, const string& caseFolder)
{
    cuttingInfo_.readDict(line3DDict.subDict(line3DObjectKeys::SURFACE_INTERSECTION_KEY), caseFolder);
}

[[nodiscard]] ItemDataSetProvider* SurfaceIntersectionLine3DSourceInfo::createDataSetProvider(const std::string& name) const
{
    return new SurfaceIntersectionLine3DDataSetProvider(
        name,
        cuttingInfo_.createCuttingSurfaceProvider()->getImplicitFunction()
    );
}

void SurfaceIntersectionLine3DSourceInfo::addItemRequirements(ItemRequirements &itemRequirements) const
{
    itemRequirements.setNeedsGhostCellsTrue();
}

void SurfaceIntersectionLine3DSourceInfo::computeAndAddToHash(size_t& hash) const
{
    cuttingInfo_.computeAndAddToHash(hash);
}

bool SurfaceIntersectionLine3DSourceInfo::operator==(const Line3DSourceInfo& other) const
{
    const auto* casted = dynamic_cast<const SurfaceIntersectionLine3DSourceInfo*>(&other);
    if (casted == nullptr)
    {
        return false;
    }
    return cuttingInfo_ == casted->cuttingInfo_;
}

} // End namespace Foam

// ************************************************************************* //
