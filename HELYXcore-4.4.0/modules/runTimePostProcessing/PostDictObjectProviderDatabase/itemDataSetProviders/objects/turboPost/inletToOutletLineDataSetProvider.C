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

#include "inletToOutletLineDataSetProvider.H"

#include "dataStructs/objects/turboPost/inletToOutletLineObjectData.H"

#include "engysTurboLines.h"

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

InletToOutletLineDataSetProvider::InletToOutletLineDataSetProvider(
    const std::string &name,
    const InletToOutletLineObjectData &dictData
)
    : ItemDataSetProvider(name)
{
    this->inletToOutletLine_ = vtkSmartPointer<engysTurboLines>::New();
    inletToOutletLine_->InletToOutletLineOn();
    inletToOutletLine_->SetLineValue(dictData.spanValue);
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void InletToOutletLineDataSetProvider::update(scalar currentTime)
{
    inletToOutletLine_->SetInputData(sources_.at(0)->getDataObjectOutput());
    inletToOutletLine_->Update();
    output_ = inletToOutletLine_->GetOutput();
}

} // End namespace Foam


// ************************************************************************* //
