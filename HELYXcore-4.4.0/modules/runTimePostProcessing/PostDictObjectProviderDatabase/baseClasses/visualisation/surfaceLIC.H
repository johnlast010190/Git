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
    (c) 2019 OpenCFD Ltd.
    (c) 2022-2023 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::SurfaceLIC

Description
    Struct that holds the Surface LIC visualisation information about an object.

SourceFiles
    surfaceLIC.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "colourMaps/colours.H"

#include "meshes/primitiveShapes/point/point.H"
#include "baseClasses/foamField.H"
#include "types/representation.H"
#include "types/backfaceStyling.H"

namespace Foam::functionObjects::runTimeVis
{
    struct SurfaceLIC
    {
        SurfaceLIC() = default;
        SurfaceLIC(const SurfaceLIC& copy) = default;

        bool visible = false;
        scalar stepSize = 1.0;
        label numberOfSteps = 10;
        word vectorField;

        void readFromDict(const dictionary& dict, const SurfaceLIC& defaultValues);

        [[nodiscard]] bool isVisible() const;

        SurfaceLIC& operator=(const SurfaceLIC&);
    };


} // End namespace Foam
