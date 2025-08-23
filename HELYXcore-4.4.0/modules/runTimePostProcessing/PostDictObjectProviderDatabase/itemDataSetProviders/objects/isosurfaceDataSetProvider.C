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
    (c) 2024-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "isosurfaceDataSetProvider.H"

#include "Utils/Utils.H"

#include "vtkContourFilter.h"


namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

IsoSurfaceDataSetProvider::IsoSurfaceDataSetProvider(const std::string& name, const IsoSurfaceObjectData& dictData)
        : ItemDataSetProvider(name)
{
    contourFilter_ = vtkSmartPointer<vtkContourFilter>::New();
    contourFilter_->ComputeScalarsOn();
    contourFilter_->ComputeNormalsOff();
    contourFilter_->GenerateTrianglesOn();
    contourFilter_->ComputeGradientsOff();

    contourFilter_->SetNumberOfContours(dictData.values.size());
    forAll(dictData.values, valuei)
    {
        contourFilter_->SetValue(valuei, dictData.values[valuei]);
    }

    contourFilter_->SetInputArrayToProcess
    (
        0,  // index: scalars(0)
        0,  // port
        0,  // connection
        vtkDataObject::FIELD_ASSOCIATION_POINTS,
        dictData.scalarField.c_str()
    );
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void IsoSurfaceDataSetProvider::update(scalar vtkNotUsed(currentTime))
{
    contourFilter_->SetInputData(sources_.at(0)->getDataSetOutput());
    contourFilter_->Update();
    output_ = contourFilter_->GetOutput();
}

} // End namespace


// ************************************************************************* //
