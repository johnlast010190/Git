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
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::PvdFileInfo

Description
    Information about a PVD geometry object.
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/itemInfo.H"
#include "dataStructs/fileSource/pvdFileData.H"
#include "itemDataSetProviders/fileSource/pvdFileDataSetProvider.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class PvdFileInfo : public ItemInfo
{
protected:
    PvdFileData data_;

public:
    explicit PvdFileInfo
    (
        const dictionary& itemDict,
        const dictionary& objectVisDict,
        const string &region,
        const string& caseFolder
    )
    :
        ItemInfo(itemDict, ItemType::Value::PVD_FILE_SOURCE, Visualisation::indexDefaults()),
        data_("", caseFolder)
    {
        id_.region = region;
        string defaultFilePath = id_.name + ".pvd";
        data_.pvdFilePath = caseFolder / itemDict.lookupOrDefault(surfaceKeys::FILE_KEY, defaultFilePath);
        readVisualisationDict(objectVisDict);
    };

    ~PvdFileInfo() override = default;

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return new PvdReaderDataSetProvider(id_.name, data_);
    };

    [[nodiscard]] PvdFileInfo* copy() const override {
        return new PvdFileInfo(*this);
    }

protected:
    [[nodiscard]] size_t computeAndReturnSubclassItemInfoHash() const override
    {
        size_t hash = ItemType::Value::PVD_FILE_SOURCE;
        data_.computeAndAddToHash(hash);
        return hash;
    }

    bool isSubclassDataEqualTo(const ItemInfo* other) const override
    {
        auto casted = dynamic_cast<const PvdFileInfo*>(other);
        if (casted == nullptr)
        {
            return false;
        }
        return data_ == casted->data_;
    };

    [[nodiscard]] bool needsGhostCellsWhenTransparent() const override {return false;};

    [[nodiscard]] bool needsGhostCellsWhenProfile() const override {return false;};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
