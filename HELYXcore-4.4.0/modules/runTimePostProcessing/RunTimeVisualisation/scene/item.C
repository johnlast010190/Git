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

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include <memory>

#include "item.H"

// VTK includes

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Item::Item
(
    const ItemInfo* info,
    int index
)
:
    info_(info),
    index(index)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Item::setActorDataSet(vtkDataSet *dataset)
{
    if (actor_ != nullptr)
    {
        actor_->processInputData(dataset);
    }
}

void Item::setActorOutline(vtkPolyData *outline)
{
    if (actor_ != nullptr)
    {
        actor_->setProcessedInputData(outline);
    }
}

void Item::registerActorDataSetToCompositer(Compositer* compositer)
{

    if (actor_ != nullptr)
    {
        compositer->registerActorPolyData(*actor_, this->info_->getId().name.c_str());
    }
}

void Item::redistributeActorDataSetIfNecessary(Compositer* compositer)
{
    if (actor_ != nullptr)
    {
        compositer->redistributeActorPolyData(*actor_);
    }
}

rtppActor* Item::getConfiguredActor()
{
    if (actor_ == nullptr)
    {
        actor_ = std::make_shared<rtppActor>(info_->getVisualisation(), index);
    }
    return actor_.get();
}

bool Item::isVisible() const
{
    return info_->isVisible();
}


} // End namespace Foam

// ************************************************************************* //
