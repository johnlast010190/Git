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
    Foam::functionObjects::runTimeVis::ProcessorBoundaryDataSetProvider

Description
    Provider for the processor boundary patches, which computes their vtk dataset

SourceFiles
    processorBoundaryDataSetProvider.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "itemDataSetProviders/mesh/patchDataSetProvider.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                   Class ProcessorBoundaryDataSetProvider Declaration
\*---------------------------------------------------------------------------*/

class ProcessorBoundaryDataSetProvider : public PatchDataSetProvider
{

public:
        ProcessorBoundaryDataSetProvider(const string& name, const string& region, const MeshAndFields& meshAndFields)
            : PatchDataSetProvider(name, region, meshAndFields) {};
        ~ProcessorBoundaryDataSetProvider() override = default;

        void update(scalar currentTime) override;

        bool acceptsProcessBoundarySources() const override {return false;};

        bool isProcessorBoundary() override {return true;};
        label getTargetProcessorBoundary() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
