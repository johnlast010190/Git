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
    Foam::functionObjects::runTimeVis::PatchDataSetProvider

Description
    Provider for the patches, which computes their vtk dataset

SourceFiles
    patchDataSetProvider.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "itemDataSetProviders/mesh/meshDataSetProvider.H"
#include "primitives/ints/lists/labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class vtkMultiPieceDataSet;

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                   Class CellZoneDataSetProvider Declaration
\*---------------------------------------------------------------------------*/

class PatchDataSetProvider : public MeshDataSetProvider
{
protected:

    bool isEmptyPatch_ = false;
    bool isCoupledPatch_ = false;
    bool isNCCPatch_ = false;

    vtkSmartPointer<vtkFloatArray> getConvertedFieldForPatch
        (
            const foamField &foamFieldName,
            const label &patchIndex,
            bool& isPointArray
        );

    vtkSmartPointer<vtkFloatArray> getConvertedFieldForEmptyPatch
        (
                const foamField &foamFieldName,
                const fvPatch &patch,
                bool& isPointArray
        );

    void addArrayToConvertedMesh
        (
            const foamField &fieldName,
            const fvPatch& patch,
            const label &patchIndex,
            vtkDataSet* convertedMesh
        );

    void update(scalar currentTime) override;

    bool acceptsProcessBoundarySources() const override
    { return true; };

public:
    PatchDataSetProvider(const string &name, const string &region, const MeshAndFields& meshAndFields);

    ~PatchDataSetProvider() override;

    bool isEmptyPatch() override
    { return isEmptyPatch_; };

    bool isCoupledPatch() override
    { return isCoupledPatch_; };

    bool isNCCPatch() override
    { return isNCCPatch_; };
private:

    void handleConditionalPatch(const fvPatch &conditionalPatch);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
