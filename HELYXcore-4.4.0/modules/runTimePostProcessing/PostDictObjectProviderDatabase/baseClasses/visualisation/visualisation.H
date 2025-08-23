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
    (c) 2019 OpenCFD Ltd.
    (c) 2022-2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::visualisation

Description
    Struct that holds visualisation information about an object.

SourceFiles
    visualisation.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "colourMaps/colours.H"

#include "meshes/primitiveShapes/point/point.H"
#include "baseClasses/foamField.H"
#include "baseClasses/foamFields.H"
#include "types/representation.H"
#include "types/backfaceStyling.H"
#include "baseClasses/visualisation/surfaceLIC.H"
#include <unordered_set>

namespace Foam::functionObjects::runTimeVis
{
    struct Visualisation
    {
        Visualisation(const dictionary& dict, const Visualisation* defaultValues);
        explicit Visualisation(const Visualisation* defaultValues);
        Visualisation(const Visualisation& copy) = default;

        static const Visualisation* solidColorDefaults();
        static const Visualisation* indexDefaults();

        Representation representation;
        bool visible = false;
        bool showActiveGIBBoundary = false;
        scalar opacity = 1.0;
        foamField colourField;
        point colour;
        point activeColour;
        point lineColour;
        scalar lineThickness = 1.0;
        dictionary supplementaryInfo;
        BackfaceStyling backfaceStyling;
        SurfaceLIC surfaceLic;

        bool showNormals = false;
        scalar normalsLength = 0.01;;
        scalar normalsRatio = 1;
        scalar normalsOpacity = 1;
        point normalsColor;

        void readFromDict(const dictionary& dict);

        [[nodiscard]] bool isVisible() const;

        [[nodiscard]] bool isTransparent() const;

        [[nodiscard]] FoamFields getRequiredFields() const;

    private:
        Visualisation() = default;

        Visualisation& operator=(const Visualisation&);

        static Visualisation createCommonDefaultValues();
        static Visualisation createSolidColorDefaultValues();
        static Visualisation createIndexDefaultValues();

        const Visualisation* const defaultVisualisation = nullptr;
    };
} // End namespace Foam
