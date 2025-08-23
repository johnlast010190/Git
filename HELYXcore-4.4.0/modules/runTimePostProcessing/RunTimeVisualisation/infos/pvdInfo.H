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
    (c) 2024-2025 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::PvdInfo

Description
    A struct with basic info about the pvd export

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "types/pvdFieldsType.H"

#include "baseClasses/foamFields.H"
#include "baseClasses/itemRequirements.H"
#include "types/compressionLevelType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{
struct SceneInfo;

struct PvdInfo {
    PvdFieldsType pvdFieldsType;
    ItemRequirements pvdRequirements;
    bool exportPvd = false;
    CompressionLevelType compressionLevel = CompressionLevelType(CompressionLevelType::Value::BALANCED);
    bool writePointData = true;
    bool oneFilePerProcess = false;
    bool removeOldObsoleteFiles = false;

    PvdInfo() = default;

    explicit PvdInfo(const dictionary& postDict);

    [[nodiscard]] const ItemRequirements& getItemRequirements() const;

    [[nodiscard]] FoamFields getAllowedFieldsForPVD(const SceneInfo& sceneInfo) const;

private:
    static const FoamFields& allFields();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
