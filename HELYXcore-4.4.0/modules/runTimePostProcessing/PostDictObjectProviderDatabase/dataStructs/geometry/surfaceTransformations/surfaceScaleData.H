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

Class
    Foam::functionObjects::runTimeVis::SurfaceScaleData

Description
    Contains the dict data necessary to scale a surface

SourceFiles
    surfaceScaleData.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "surfaceTransformationData.H"
#include "primitives/Vector/vector/vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations

namespace Foam::functionObjects::runTimeVis
{

struct SurfaceScaleData : public SurfaceTransformationData {
    Vector<scalar> origin;
    Vector<scalar> scale;

    explicit SurfaceScaleData(const dictionary& transformationDict)
    {
        SurfaceScaleData::readDict(transformationDict);
    }

    void readDict(const dictionary& transformationDict) override;

    bool operator==(const SurfaceTransformationData& other) const override;

    void computeAndAddToHash(size_t& hash) const override;

    [[nodiscard]] std::unique_ptr<SurfaceTransformationData> copy() const override;

    void addToVtkTransform(vtkTransform*) const override;

private:
    SurfaceScaleData() = default;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
