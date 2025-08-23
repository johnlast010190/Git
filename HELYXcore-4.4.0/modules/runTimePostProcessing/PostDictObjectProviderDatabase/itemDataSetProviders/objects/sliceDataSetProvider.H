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
    Foam::functionObjects::runTimeVis::SliceDataSetProvider

Description
    Provider for the slice 3d objects, which computes their vtk dataset

SourceFiles
    SliceDataSetProvider.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "dataStructs/objects/sliceObjectData.H"
#include "itemDataSetProviders/itemDataSetProvider.H"
#include "itemDataSetProviders/objects/cutting/cuttingSurfaceProviders.H"
#include "hash/hasher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class engysSliceDataSet;


namespace Foam::functionObjects::runTimeVis
{

class ReferenceFrames;
class ReferenceFrame;

/*---------------------------------------------------------------------------*\
                   Class SliceDataSetProvider Declaration
\*---------------------------------------------------------------------------*/

class SliceDataSetProvider : public ItemDataSetProvider
{
private:
    const std::shared_ptr<const CuttingSurfaceProvider> cutProvider_;
    vtkSmartPointer<engysSliceDataSet> sliceAlgorithm_;

    vtkSmartPointer<vtkImplicitFunction> implicitFunction_;
    const ReferenceFrame* referenceFrame_;

public:
    SliceDataSetProvider
    (
            const std::string& name,
            const SliceObjectData& dictData,
            const std::shared_ptr<const CuttingSurfaceProvider>& cutProvider,
            const ReferenceFrames& referenceFrames
    );
    ~SliceDataSetProvider() override = default;


private:
    void update(scalar currentTime) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
