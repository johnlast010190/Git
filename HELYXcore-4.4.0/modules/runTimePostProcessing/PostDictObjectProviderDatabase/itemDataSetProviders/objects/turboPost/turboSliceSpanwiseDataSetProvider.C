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

\*---------------------------------------------------------------------------*/

#include <vtkAppendPolyData.h>
#include "turboSliceSpanwiseDataSetProvider.H"

#include "dataStructs/objects/turboPost/turboSliceSpanwiseObjectData.H"

#include "engysParametricSlicer.h"
#include "engysParametricBackgroundMesh.h"

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

TurboSliceSpanwiseDataSetProvider::TurboSliceSpanwiseDataSetProvider(
    const std::string &name,
    const TurboSliceSpanwiseObjectData &dictData
)
    : ItemDataSetProvider(name)
{
    this->turboSlicer_ = vtkSmartPointer<engysParametricSlicer>::New();
    turboSlicer_->SetUnravel(dictData.unravel);
    turboSlicer_->SetValueToSlice(dictData.streamSlice);
    turboSlicer_->SetArrayToSlice(engysParametricBackgroundMesh::GetStreamArrayName());
    turboSlicer_->SetArrayToBound(engysParametricBackgroundMesh::GetSpanArrayName());
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void TurboSliceSpanwiseDataSetProvider::update(scalar currentTime)
{
    vtkDataObject *multiBlock = sources_.at(0)->getDataObjectOutput();
    turboSlicer_->SetInputData(multiBlock);
    turboSlicer_->Update();
    output_ = turboSlicer_->GetOutput();
}

} // End namespace Foam


// ************************************************************************* //
