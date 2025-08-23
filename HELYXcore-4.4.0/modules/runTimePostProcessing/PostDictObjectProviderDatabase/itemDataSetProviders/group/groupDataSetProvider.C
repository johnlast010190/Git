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
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "groupDataSetProvider.H"

#include "engysAppendFilter.h"

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

GroupDataSetProvider::GroupDataSetProvider(const std::string& name, const GroupData& vtkNotUsed(dictData))
        : ItemDataSetProvider(name)
{
    append_ = vtkSmartPointer<engysAppendFilter>::New();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void GroupDataSetProvider::update(scalar vtkNotUsed(currentTime))
{
    if (anySourceHasGhostCells())
    {
        allocateSourcesGhostCells();
    }

    append_->RemoveAllInputs();
    for (ItemDataSetProvider* source : sources_)
    {
        append_->AddInputData(source->getDataSetOutput());
    }
    append_->Update();
    output_ = append_->GetOutput();
}

bool GroupDataSetProvider::anySourceHasGhostCells()
{
    bool needsGhosts = false;
    for (ItemDataSetProvider* source : sources_)
    {
        if (source->getDataSetOutput()->HasAnyGhostCells())
        {
            needsGhosts = true;
            break;
        }
    }
    return needsGhosts;
}

void GroupDataSetProvider::allocateSourcesGhostCells()
{
    for (ItemDataSetProvider* source : sources_)
    {
        if (!source->getDataSetOutput()->GetCellGhostArray())
        {
            source->getDataSetOutput()->AllocateCellGhostArray();
        }
    }
}

} // End namespace


// ************************************************************************* //
