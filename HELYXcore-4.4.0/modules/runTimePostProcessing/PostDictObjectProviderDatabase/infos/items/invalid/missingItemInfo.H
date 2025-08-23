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
    (c) 2023-2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::MissingPostDictItemInfo

Description
    Information about items that are part of a scene but don't exist in the base model.

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/invalid/invalidItemInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class MissingPostDictItemInfo : public InvalidItemInfo
{
public:
    explicit MissingPostDictItemInfo(const Id& id)
    :
        InvalidItemInfo(id)
    {};
    ~MissingPostDictItemInfo() override = default;

    [[nodiscard]] word getInvalidReason() const override {return std::string("was not defined in the postDict or runtimeVisualisationDict");};

    [[nodiscard]] MissingPostDictItemInfo* copy() const override {
        return new MissingPostDictItemInfo(*this);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
