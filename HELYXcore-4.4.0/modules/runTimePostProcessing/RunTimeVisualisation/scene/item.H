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
    Foam::functionObjects::runTimeVis::Item

Description
    A class that combines the item info, dataset provider, and actor.

SourceFiles
    item.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "rendering/compositers/compositer.H"
#include "infos/items/itemInfo.H"
#include "rendering/rtppActor.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                          Class Item Declaration
\*---------------------------------------------------------------------------*/

class Item
{
private:


protected:

        const ItemInfo* info_;
        std::shared_ptr<rtppActor> actor_ = nullptr;
        int index = 0;

public:

    // Constructors

        explicit Item
        (
            const ItemInfo* info,
            int index
        );

    // Member Functions

        void setActorDataSet(vtkDataSet *dataset);
        void setActorOutline(vtkPolyData *outline);

        void registerActorDataSetToCompositer(Compositer* compositer);
        void redistributeActorDataSetIfNecessary(Compositer* compositer);

        rtppActor* getConfiguredActor();

        const ItemInfo* getInfo() {return info_;};

        const foamField& getColorField() {return info_->getColorField();};

        [[nodiscard]] bool isVisible() const;

        //- No copy assignment
        Item& operator=(const Item&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
