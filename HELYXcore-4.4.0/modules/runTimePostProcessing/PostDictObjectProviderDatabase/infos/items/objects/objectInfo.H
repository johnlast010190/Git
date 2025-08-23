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
    Foam::functionObjects::runTimeVis::ObjectInfo

Description
    A base class for all renderable runTimeVisualisation object infos.
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/itemInfo.H"
#include "postDict/postDictKeys.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                          Class surface Declaration
\*---------------------------------------------------------------------------*/

class ObjectInfo : public ItemInfo
{
private:


protected:

    Id source_;

    ObjectInfo() : ItemInfo(Visualisation::solidColorDefaults()), source_() {};

public:

    // Constructors

    explicit ObjectInfo
    (
        const dictionary& dict
    ) : ItemInfo(dict, Visualisation::solidColorDefaults()),
    source_(dict.subDict(objectKeys::SOURCE_KEY))
    {};
    ~ObjectInfo() override = 0;

    [[nodiscard]] std::vector<Id> getSources() const override
    {
        std::vector<Id> sources;
        sources.push_back(source_);
        return sources;
    };

    [[nodiscard]] bool isObject() const override {return true;};

protected:
    template<class T, class Data>
    inline bool isTypeAndDataEqualTo(const ItemInfo* other, const Data& data) const
    {
        const T* casted = dynamic_cast<const T*>(other);
        if (casted == nullptr)
        {
            return false;
        }
        return data == casted->getData();
    };
};
ObjectInfo::~ObjectInfo() = default;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
