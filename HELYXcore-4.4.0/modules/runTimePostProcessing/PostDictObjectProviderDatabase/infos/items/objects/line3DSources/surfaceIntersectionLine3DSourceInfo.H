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
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::SurfaceIntersectionLine3DSourceData

Description
    Contains the info for a 3D line from an intersection with another surface

SourceFiles
    surfaceIntersectionLine3DSourceInfo.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/objects/line3DSources/line3DSourceInfo.H"
#include "db/dictionary/dictionary.H"
#include "infos/items/objects/cutting/cuttingInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class SurfaceIntersectionLine3DSourceInfo : public Line3DSourceInfo
{
public:

    explicit SurfaceIntersectionLine3DSourceInfo(const dictionary& line3DDict, const string& caseFolder);

    ~SurfaceIntersectionLine3DSourceInfo() override = default;

    void readDict(const dictionary& line3DDict, const string& caseFolder) override;

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider(const std::string& name) const override;

    void addItemRequirements(ItemRequirements& itemRequirements) const override;

    void computeAndAddToHash(size_t& hash) const override;

    bool operator==(const Line3DSourceInfo& other) const override;
private:
    CuttingInfo cuttingInfo_;
};

} // End namespace Foam

// ************************************************************************* //
