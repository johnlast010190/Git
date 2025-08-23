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

\*---------------------------------------------------------------------------*/

#include "importDataDataSetProvider.H"

#include "vtkMultiProcessController.h"

namespace Foam::functionObjects::runTimeVis
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ImportDataDataSetProvider::ImportDataDataSetProvider(const std::string& name, const ImportDataData& dictData)
: ExternalItemDataSetProvider(name)
{
    reader_ = vtkSmartPointer<engysExternalDataSetReader>::New();
    reader_->SetFilePath(dictData.path.c_str());
    reader_->SetTimeValue(dictData.timeStep);
    reader_->SetFileFormat(dictData.format);
    reader_->SetController(vtkMultiProcessController::GetGlobalController());
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void ImportDataDataSetProvider::update(scalar vtkNotUsed(currentTime))
{
    if (!output_.Get())
    {
        reader_->Update();
        output_ = reader_->GetOutput();
    }
}

} // End namespace


// ************************************************************************* //
