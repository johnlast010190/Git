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
    (c) 2019 OpenCFD Ltd.
    (c) 2025 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::ItemRequirements

Description
    Class that stores the requirements from downstream, so that the upstream
    objects can fulfill them

SourceFiles
    <itemRequirements.C>

\*---------------------------------------------------------------------------*/
#pragma once

#include "itemRequirements.H"

#include "foamFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class ItemRequirements
{
public:

    // Constructors
    ItemRequirements() = default;

    // Downstream items (items closer to what is shown in the scene) call this in each upstream
    // item (items closer to the mesh) so that they can fulfill their needs
    void mergeFromDownstream(const ItemRequirements& downstreamRequirements);

    const FoamFields& getRequiredFields() const
    {return requiredFields_; }
    FoamFields& getWritableRequiredFields()
    {return requiredFields_; }

    void addToRequiredFields(const foamField& field);

    void setNeedsGhostCellsTrue()
    { this->needsGhostCells_ = true; }
    bool getNeedsGhostCells() const
    { return this->needsGhostCells_; }

    ~ItemRequirements() = default;

private:
    FoamFields requiredFields_;

    bool needsGhostCells_ = false;
};

} // End namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
