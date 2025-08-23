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
    (c) 2022-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "pvdFileDataSetProvider.H"

#include "vtkMultiProcessController.h"
#include "engysPVDReader.h"
#include "vtkUnstructuredGrid.h"
#include "vtkImageData.h"

namespace Foam::functionObjects::runTimeVis
{

PvdReaderDataSetProvider::PvdReaderDataSetProvider(const std::string &name, PvdFileData data)
    : ItemDataSetProvider(name), data_(std::move(data))
{
    if (data_.pvdFilePath.empty())
    {
        FatalErrorInFunction << "No pvd file specified for " << name_ << abort(FatalError);
    }

    reader_ = vtkSmartPointer<engysPVDReader>::New();
    reader_->SetFileName(data_.pvdFilePath.c_str());
    reader_->SetController(vtkMultiProcessController::GetGlobalController());
}

void PvdReaderDataSetProvider::update(scalar currentTime)
{
    if (reader_->FileExists())
    {
        reader_->Modified(); // Necessary so that it reads the available timesteps in the pvd again
        reader_->UpdateTimeStep(currentTime, -1, 0, 0);
        output_ = reader_->GetOutputAsDataSet();
    }
    else
    {
        FatalErrorInFunction << "Pvd file " << data_.pvdFilePath << " was not found!" << abort(FatalError);
        output_ = nullptr;
    }

    if (output_ == nullptr)
    {
        output_ = vtkSmartPointer<vtkUnstructuredGrid>::New();
    }
}

} // End namespace Foam


// ************************************************************************* //
