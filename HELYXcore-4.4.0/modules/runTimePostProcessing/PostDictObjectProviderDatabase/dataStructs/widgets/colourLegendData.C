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

// OpenFOAM includes
#include "colourLegendData.H"

namespace Foam::functionObjects::runTimeVis
{

ColourLegendData::ColourLegendData() :
    title("default"),
    visible(true),
    labelFormat("%6.2f"),
    numberOfLabels(5),
    location(),
    coordinates(0.0, 0.0),
    thickness(16),
    length(0.5),
    vertical(true),
    autoLabels(true),
    showTicks(false)
    {}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ColourLegendData::ColourLegendData
(
    const word& name,
    const dictionary& legendDict
) : ColourLegendData(name, legendDict, ColourLegendData(), true)
{
}

ColourLegendData::ColourLegendData
(
    const word& name,
    const dictionary& legendDict,
    const ColourLegendData& defaultColorLegendData,
    bool allowNoTitle
) :
name(name),
coordinates(0,0)
{
    if (allowNoTitle) {
        title = legendDict.lookupOrDefault<string>(colourLegendsKeys::TITLE_KEY, defaultColorLegendData.title);
    } else {
        title = legendDict.lookup(colourLegendsKeys::TITLE_KEY);
    }
    visible = legendDict.lookupOrDefault<bool>(visualisationKeys::VISIBLE_KEY, defaultColorLegendData.visible);
    thickness = legendDict.lookupOrDefault<label>(colourLegendsKeys::THICKNESS_KEY, defaultColorLegendData.thickness);
    length = legendDict.lookupOrDefault<scalar>(colourLegendsKeys::HEIGHT_KEY, defaultColorLegendData.length);
    location = legendDict.lookupOrDefault<Location>(colourLegendsKeys::LOCATION_KEY, defaultColorLegendData.location);
    if (location.IsCustom())
    {
        coordinates = legendDict.lookupOrDefault<Tuple2<scalar, scalar>>(colourLegendsKeys::COORDINATES_KEY, defaultColorLegendData.coordinates);
    }

    vertical = legendDict.lookupOrDefault<bool>(colourLegendsKeys::VERTICAL_KEY, defaultColorLegendData.vertical);
    labelFormat = legendDict.lookupOrDefault<string>(colourLegendsKeys::LABEL_FORMAT_KEY, defaultColorLegendData.labelFormat);
    numberOfLabels = legendDict.lookupOrDefault<label>(colourLegendsKeys::NUMBER_OF_LABELS_KEY, defaultColorLegendData.numberOfLabels);
    autoLabels = legendDict.lookupOrDefault<bool>(colourLegendsKeys::AUTOMATIC_LABELS_KEY, defaultColorLegendData.autoLabels);
    showTicks = legendDict.lookupOrDefault<bool>(colourLegendsKeys::SHOW_TICKS_KEY, defaultColorLegendData.showTicks);

    if (legendDict.found(colourLegendsKeys::FONT_DICT_KEY))
    {
        const dictionary& labelFontDict = &legendDict.subDict(colourLegendsKeys::FONT_DICT_KEY);
        font.readDict(labelFontDict, defaultColorLegendData.font);
    }
    else
    {
        font = defaultColorLegendData.font;
    }
}

} // End namespace


// ************************************************************************* //
