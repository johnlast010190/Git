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
    Foam::functionObjects::runTimeVis::SurfaceDataSetProvider

Description
    Base class for the dataset providers of the surfaces

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "itemDataSetProviders/itemDataSetProvider.H"
#include "dataStructs/geometry/surfaceTransformationsData.H"
#include "Utils/ParallelUtils.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
               Class SurfaceDataSetProvider Declaration
\*---------------------------------------------------------------------------*/

class SurfaceDataSetProvider : public ItemDataSetProvider
{
protected:

    SurfaceDataSetProvider(const std::string &name, const SurfaceTransformationsData& transformationsData);
    explicit SurfaceDataSetProvider(const std::string &name);

    bool initialised_ = false;

    virtual void initialiseSurface() = 0;

    void update(scalar currentTime) override;

    [[nodiscard]] vtkSmartPointer<vtkPolyData> getDataSetOutline() override;

    ~SurfaceDataSetProvider() override = 0;

private:
    SurfaceTransformationsData transformationsData_;
    double originalBounds_[6];
    vtkSmartPointer<vtkPolyData> outline_;
    void transformSurface();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
