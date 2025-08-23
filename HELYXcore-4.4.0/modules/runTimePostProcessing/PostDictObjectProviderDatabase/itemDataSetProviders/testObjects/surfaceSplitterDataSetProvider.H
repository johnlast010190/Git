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

Class
    Foam::functionObjects::runTimeVis::SurfaceSplitterDataSetProvider

Description
    Provider for the clip 3d objects, which computes their vtk dataset

SourceFiles
    surfaceSplitterDataSetProvider.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "dataStructs/testObjects/surfaceSplitterObjectData.h"
#include "itemDataSetProviders/itemDataSetProvider.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
              Class SurfaceSplitterDataSetProvider Declaration
\*---------------------------------------------------------------------------*/

class SurfaceSplitterDataSetProvider : public ItemDataSetProvider
{
public:
    SurfaceSplitterDataSetProvider
    (
            const std::string& name,
            const SurfaceSplitterObjectData& dictData
    );
    ~SurfaceSplitterDataSetProvider() override = default;

private:
    SurfaceSplitterObjectData data_;

    void update(scalar currentTime) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
