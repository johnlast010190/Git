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
    Foam::functionObjects::runTimeVis::SurfaceTransformationData

Description
    Contains the dict data necessary to transform a surface with one transformation step

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "db/dictionary/dictionary.H"
#include "vtkSmartPointer.h"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class vtkTransform;

namespace Foam::functionObjects::runTimeVis
{

struct SurfaceTransformationData {

    virtual void readDict(const dictionary& transformationDict) = 0;

    virtual bool operator==(const SurfaceTransformationData& other) const = 0;

    virtual void computeAndAddToHash(size_t& hash) const = 0;

    virtual void addToVtkTransform(vtkTransform*) const = 0;

    [[nodiscard]] virtual std::unique_ptr<SurfaceTransformationData> copy() const = 0;

    virtual ~SurfaceTransformationData() = default;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
