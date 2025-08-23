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
    Foam::functionObjects::runTimeVis::WidgetsInfo

Description
    Information about the widgets

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "dataStructs/widgets/axisWidgetData.H"
#include "dataStructs/widgets/colourLegendsData.H"
#include "dataStructs/widgets/gridWidgetData.H"
#include "dataStructs/widgets/logoWidgetData.H"
#include "dataStructs/widgets/timestepWidgetData.H"
#include "dataStructs/widgets/vectorWidgetData.H"
#include "dataStructs/widgets/referenceFramesWidgetData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct WidgetsInfo
{
protected:
    const dictionary& widgetsSubDict_;

public:
    AxisWidgetData axisWidgetData;
    ColourLegendsData colourLegendsData;
    GridWidgetData gridWidgetData;
    LogoWidgetData logoWidgetData;
    TimestepWidgetData timestepWidgetData;
    VectorWidgetData vectorWidgetData;
    ReferenceFramesWidgetData referenceFramesWidgetData;


    // Construct from dictionary
    explicit WidgetsInfo(const dictionary& sceneDict) :
        widgetsSubDict_(GET_OPTIONAL_DICTIONARY(sceneDict, sceneKeys::WIDGETS_DICT_KEY)),
        axisWidgetData(widgetsSubDict_),
        colourLegendsData(sceneDict),
        gridWidgetData(widgetsSubDict_),
        logoWidgetData(widgetsSubDict_),
        timestepWidgetData(widgetsSubDict_),
        vectorWidgetData(widgetsSubDict_),
        referenceFramesWidgetData(sceneDict)
    {};
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
