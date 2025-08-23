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
    Foam::functionObjects::runTimeVis::FieldSamplingFunctionObjectInfo

Description
    Information a line geometry object.
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/itemInfo.H"
#include "dataStructs/fileSource/fieldSamplingPvdFileData.H"
#include "itemDataSetProviders/fileSource/pvdFileDataSetProvider.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class FieldSamplingFunctionObjectInfo : public ItemInfo
{
protected:
    FieldSamplingPvdFileData data_;

public:
    explicit FieldSamplingFunctionObjectInfo
    (
        const string& fileSourceItemName,
        const word& region,
        const string& caseFolder,
        const dictionary& functionObjectDict
    )
    :
        ItemInfo(Visualisation::solidColorDefaults()),
        data_(functionObjectDict, caseFolder)
    {
        id_.name = fileSourceItemName;
        id_.region = region;
        id_.type = ItemType(ItemType::FIELD_SAMPLING_FILE_SOURCE);
    };

    ~FieldSamplingFunctionObjectInfo() override = default;

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return new PvdReaderDataSetProvider(id_.name, data_);
    };

    [[nodiscard]] FieldSamplingFunctionObjectInfo* copy() const override {
        return new FieldSamplingFunctionObjectInfo(*this);
    }

protected:
    [[nodiscard]] size_t computeAndReturnSubclassItemInfoHash() const override
    {
        size_t hash = 0;
        data_.computeAndAddToHash(hash);
        return hash;
    }

    bool isSubclassDataEqualTo(const ItemInfo* other) const override
    {
        auto casted = dynamic_cast<const FieldSamplingFunctionObjectInfo*>(other);
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

} // End namespace runTimeVis
// End namespace functionObjects
// End namespace Foam

// ************************************************************************* //
