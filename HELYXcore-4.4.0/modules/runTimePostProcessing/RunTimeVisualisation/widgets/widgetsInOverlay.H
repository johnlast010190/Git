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
    Foam::functionObjects::runTimeVis::WidgetsInOverlay

Description
    A class that holds all widgets whose actors are meant to be overlaid over the 3D actors

SourceFiles
    widgetsInOverlay.C

\*---------------------------------------------------------------------------*/

#pragma once

#include <memory>
#include "primitives/Scalar/scalar/scalar.H"
#include "db/Time/Time.H"

#include "widgets/overlay/colorLegends.H"
#include "widgets/overlay/vectorWidget.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
class vtkRenderer;

namespace Foam::functionObjects::runTimeVis
{

struct WidgetsInfo;
class AxesWidget;
class LogoWidget;
class TimestepWidget;
class ColourLookupTableProvider;

/*---------------------------------------------------------------------------*\
                        Class WidgetsInOverlay Declaration
\*---------------------------------------------------------------------------*/

class WidgetsInOverlay
{
private:
    vtkRenderer* renderer_;

    std::unique_ptr<AxesWidget> axesWidget_;
    std::unique_ptr<ColorLegends> colorLegends_;
    std::unique_ptr<LogoWidget> logoWidget_;
    std::unique_ptr<TimestepWidget> timestepWidget_;
    std::unique_ptr<VectorWidget> vectorWidget_;

public:
    explicit WidgetsInOverlay(const WidgetsInfo& widgetsInfo);

    void setRenderer(vtkRenderer* overlayRenderer) {renderer_ = overlayRenderer;};

    void addVisibleWidgetsToRenderer();

    void updateVisibleWidgets(
            scalar currentTimeValue,
            ColourLookupTableProvider& colourLutProvider,
            const ColourMaps& colourMaps,
            const FoamMeshes& meshes,
            const ExternalFields& externalFields,
            const Time& runTime
    );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
