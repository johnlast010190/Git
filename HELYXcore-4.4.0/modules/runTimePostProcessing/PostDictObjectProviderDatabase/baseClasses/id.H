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
    (c) 2019 OpenCFD Ltd.
    (c) 2022 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::id

Description
    Struct that holds basic information about an object.

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include <utility>

#include "primitives/strings/word/word.H"
#include "db/dictionary/dictionary.H"
#include "Utils/Utils.H"
#include "postDict/postDictKeys.H"
#include "types/itemType.H"

namespace Foam::functionObjects::runTimeVis
{
    struct Id
    {
        string name = "";
        word region = "";
        ItemType type;

        // Constructors

        Id() = default;

        virtual ~Id() = default;

        explicit Id(const dictionary& dict)
        :
        type(dict.lookup<word>(idKeys::TYPE_KEY))
        {
            readDictTypeless(dict);
        }

        Id(const dictionary& dict, ItemType::Value type)
        :
        type(type)
        {
            readDictTypeless(dict);
        }

        Id(string  name, const string&  region, const ItemType& type)
            : name(std::move(name)), region(region), type(type)
        {}

        Id(string  name, const string&  region, ItemType::Value type)
            : name(std::move(name)), region(region), type(type)
        {}

        void readDictTypeless(const dictionary& dict)
        {
            if (dict.found(idKeys::NAME_KEY))
            {
                name = dict.lookup<string>(idKeys::NAME_KEY);
            }
            else
            {
                name = Utils::getSubdictName(dict);
            }
            region = dict.lookupOrDefault<word>(idKeys::REGION_KEY, regionKeys::DEFAULT_REGION_KEY);
        }

        bool operator==(const Id& other) const
        {
            return
            (
                name == other.name &&
                region == other.region &&
                type == other.type
            );
        }

        bool operator!=(const Id& other) const
        {
            return
            (
                name != other.name ||
                region != other.region ||
                type != other.type
            );
        }

        friend Ostream& operator<<(Ostream& os, const Id& ID)
        {
            os << ID.name
                << " (region: " << ID.region
                << ", type: " << ID.type
                << ")";
            return os;
        };
    };


    // Hashing function for Id, taken from here:
    // https://stackoverflow.com/a/17017281
    // Allows use of Id as key in std::unordered_map
    struct IdHasher
    {
        std::size_t operator() (const Id& ID) const
        {
            // Compute individual hash values for first, second and third
            // http://stackoverflow.com/a/1646913/126995
            size_t res = 17;
            res = res * 31 + std::hash<std::string>()( ID.name );
            res = res * 31 + std::hash<int>()( ID.type.getValue() );
            res = res * 31 + std::hash<std::string>()( ID.region );
            return res;
        }
    };
} // End namespace Foam
