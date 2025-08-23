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
    (c) 2022-2025 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::ItemInfo

Description
    A base class for all renderable runTimeVisualisation item infos.
    Contains the data from the dictionary

SourceFiles
    ItemInfo.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "itemDataSetProviders/itemDataSetProvider.H"
#include "baseClasses/id.H"
#include "baseClasses/visualisation/visualisation.H"
#include "baseClasses/itemRequirements.H"

#include "hash/hasher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                          Class ItemInfo Declaration
\*---------------------------------------------------------------------------*/

class ItemInfo
{
private:


protected:

    Id id_;
    Visualisation visualisation_;

// Protected Member Functions
    [[nodiscard]] virtual size_t computeAndReturnSubclassItemInfoHash() const = 0;

    [[nodiscard]] virtual bool isSubclassDataEqualTo(const ItemInfo* other) const = 0;

    void checkObjectType(ItemType::Value expectedType) const;

    explicit ItemInfo(const Visualisation* defaultVisualisation);

public:

    // Constructors

    explicit ItemInfo
    (
        const dictionary& dict,
        const Visualisation* defaultVisualisation
    );

    ItemInfo
    (
        const dictionary& dict,
        ItemType::Value type,
        const Visualisation* defaultVisualisation
    );


    //- Destructor
    virtual ~ItemInfo();


    // Member Functions

    [[nodiscard]] virtual ItemDataSetProvider* createDataSetProvider() const = 0;

    void readVisualisationDict(const dictionary& dict);

    [[nodiscard]] virtual std::vector<Id> getSources() const {return {};};

    [[nodiscard]] const Id& getId() const {return id_;};

    [[nodiscard]] virtual bool isObject() const {return false;};
    virtual void updateObjectData(const dictionary& dict) {};

    [[nodiscard]] const Visualisation& getVisualisation() const {return visualisation_;};

    [[nodiscard]] virtual bool isExternal() const
    {
        return false;
    };

    [[nodiscard]] bool isVisible() const;
    [[nodiscard]] bool isTransparent() const;
    [[nodiscard]] const foamField& getColorField() const;

    [[nodiscard]] ItemRequirements getItemRequirements() const;

    [[nodiscard]] virtual bool isMeshBoundary() const {return false;};
    [[nodiscard]] virtual bool isProcessBoundary() const {return false;};
    [[nodiscard]] virtual bool isExternalBoundary() const {return false;};
    [[nodiscard]] virtual bool isInternalBoundary() const {return false;};

    [[nodiscard]] virtual ItemInfo* copy() const = 0;

    [[nodiscard]] virtual bool isValid() const {return true;};
    [[nodiscard]] virtual word getInvalidReason() const {return "";};


    [[nodiscard]] size_t computeAndReturnItemInfoHash() const {
        size_t hash = computeAndReturnSubclassItemInfoHash();
        hasher::hash_combine(hash, id_.type.getValue());
        return hash;
    };

    [[nodiscard]] bool isDataEqualTo(const ItemInfo* other) const;

protected:
    // True for items that need ghost cells when the item is transparent, like volume mesh, clip, threshold
    // False for items like patches, slices, glyphs
    [[nodiscard]] virtual bool needsGhostCellsWhenTransparent() const = 0;

    [[nodiscard]] virtual bool needsGhostCellsWhenProfile() const {return true;};

    virtual void addExtraItemRequirements(ItemRequirements& itemRequirements) const {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
