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
    Foam::functionObjects::runTimeVis::ItemDataSetProvider

Description
    A base class for all dataset providers, which compute the dataset for an item

SourceFiles
    ItemDataSetProvider.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "baseClasses/itemRequirements.H"
#include "baseClasses/foamField.H"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"
#include <unordered_set>

#include "vtkUnstructuredGrid.h"
#include "vtkPolyData.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations


namespace Foam::functionObjects::runTimeVis
{
class InstanceServices;

/*---------------------------------------------------------------------------*\
                   Class ItemDataSetProvider Declaration
\*---------------------------------------------------------------------------*/

class ItemDataSetProvider
{
private:
    ItemRequirements itemRequirements_;

    // TODO this can go into the global objects that should be implemented
    FoamFields nonExistentRequiredFields_;

    void addItemRequirementsToSources();

protected:

    Foam::label latestProcessedTimeIndex_;

    vtkSmartPointer<vtkDataObject> output_;

    std::vector<ItemDataSetProvider *> sources_;
    InstanceServices* instanceServices_;

    const std::string name_;

    virtual const std::vector<ItemDataSetProvider*>& getUsedSources() {return sources_; };

    void addItemRequirements(const ItemRequirements& itemRequirements);

    virtual void addExtraRequirementsDueToUpstreamRequirements(
        ItemRequirements &ownRequirements,
        const ItemRequirements &upstreamRequirements
    ) const
    {};

    [[nodiscard]] const ItemRequirements &getItemRequirements()
    { return itemRequirements_; };

    void addToNonExistingFields(const foamField& field);

    virtual void update(scalar currentTime)
    {};

public:
    explicit ItemDataSetProvider(std::string name);
    virtual ~ItemDataSetProvider() = 0;

    void setInstanceServicesForThisAndSources(InstanceServices* instanceServices);

    void addItemRequirementsForThisAndSources(const ItemRequirements& itemRequirements);

    [[nodiscard]] bool needsToUpdate(Foam::label currentTimeIndex) const
    { return latestProcessedTimeIndex_ != currentTimeIndex; };

    // updates the data if the time changed
    void updateIfNecessary(Foam::label currentTimeIndex, Foam::scalar currentTime);

    [[nodiscard]] vtkSmartPointer<vtkDataSet> getDataSetOutput(unsigned int fallbackBlock = 0);

    [[nodiscard]] vtkSmartPointer<vtkDataObject> getDataObjectOutput()
    { return output_; };

    [[nodiscard]] virtual vtkSmartPointer<vtkPolyData> getDataSetOutline();

    void addSource(ItemDataSetProvider *source)
    { sources_.push_back(source); };

    [[nodiscard]] virtual bool acceptsMeshBoundarySources() const
    { return false; };

    [[nodiscard]] virtual bool acceptsProcessBoundarySources() const
    { return false; };

    [[nodiscard]] virtual string getName()
    { return ""; };

    [[nodiscard]] virtual bool isProcessorBoundary()
    { return false; };

    [[nodiscard]] virtual bool isEmptyPatch()
    { return false; };

    [[nodiscard]] virtual bool isCoupledPatch()
    { return false; };

    [[nodiscard]] virtual bool isNCCPatch()
    { return false; };

    [[nodiscard]] virtual label getTargetProcessorBoundary()
    { return -1; };

    //- No copy construct
    ItemDataSetProvider(const ItemDataSetProvider &) = delete;

    //- No copy assignment
    void operator=(const ItemDataSetProvider &) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
