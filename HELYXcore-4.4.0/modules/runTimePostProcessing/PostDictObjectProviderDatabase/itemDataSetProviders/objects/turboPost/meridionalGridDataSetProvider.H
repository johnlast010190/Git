/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : Dev
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
    (c) 2023 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::MeridionalGridDataSetProvider

Description
    Provider for the turbo background mesh/meridional plot post-processing objects,
    which is the base for the other 4 turbo post objects

SourceFiles
    turbopropDataSetProvider.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "itemDataSetProviders/itemDataSetProvider.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class engysTurboPostBase;
class engysMeridionalPlot;

namespace Foam::functionObjects::runTimeVis
{
struct MeridionalGridObjectData;

/*---------------------------------------------------------------------------*\
               Class MeridionalGridDataSetProvider Declaration
\*---------------------------------------------------------------------------*/

class MeridionalGridDataSetProvider : public ItemDataSetProvider
{
private:

protected:

    vtkSmartPointer<engysTurboPostBase> turboBackgroundMesh_;
    vtkSmartPointer<engysMeridionalPlot> meridionalPlot_;

    int inletPatches, outletPatches, hubPatches, shroudPatches;

    void update(scalar currentTime) override;

public:
    MeridionalGridDataSetProvider(const std::string& name, const MeridionalGridObjectData& dictData);
    ~MeridionalGridDataSetProvider() override = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
