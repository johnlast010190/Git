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
    (c) 2019 OpenCFD Ltd.
    (c) 2022-2025 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::visualisation

Description
    Struct that holds visualisation information about an object.

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#include "visualisation.H"

#include "postDict/postDictKeys.H"

#include "primitives/strings/word/word.H"
#include "db/dictionary/dictionary.H"

namespace Foam::functionObjects::runTimeVis
{

Visualisation::Visualisation(const dictionary& dict, const Visualisation* defaultValues)
: defaultVisualisation(defaultValues)
{
    readFromDict(dict);
}

Visualisation::Visualisation(const Visualisation* defaultValues)
: defaultVisualisation(defaultValues)
{
    *this = (*defaultValues);
}

void Visualisation::readFromDict(const dictionary& dict)
{
    if (!defaultVisualisation)
    {
        FatalErrorInFunction
            << "A default visualisation was not set before trying to read dictionary " << dict
            << exit(FatalError);
        return;
    }

    representation = dict.lookupOrDefault<Representation>(visualisationKeys::REPRESENTATION_KEY, defaultVisualisation->representation);
    visible = dict.lookupOrDefault<bool>(visualisationKeys::VISIBLE_KEY, defaultVisualisation->visible);
    opacity = dict.lookupOrDefault<scalar>(visualisationKeys::OPACITY_KEY, defaultVisualisation->opacity);
    colourField = foamField(dict.lookupOrDefault<string>(visualisationKeys::COLOR_FIELD_KEY, defaultVisualisation->colourField));
    colour = dict.lookupOrDefault<point>(visualisationKeys::COLOR_KEY, defaultVisualisation->colour);
    lineColour = dict.lookupOrDefault<point>(visualisationKeys::LINE_COLOR_KEY, defaultVisualisation->lineColour);
    lineThickness = dict.lookupOrDefault<scalar>(visualisationKeys::LINE_THICKNESS_KEY, defaultVisualisation->lineThickness);
    backfaceStyling = dict.lookupOrDefault<BackfaceStyling>(visualisationKeys::BACKFACE_STYLING_KEY, defaultVisualisation->backfaceStyling);

    showNormals = dict.lookupOrDefault<bool>(visualisationKeys::SHOW_NORMALS_KEY, defaultVisualisation->showNormals);
    normalsLength = dict.lookupOrDefault<scalar>(visualisationKeys::NORMALS_LENGTH_KEY, defaultVisualisation->normalsLength);
    normalsRatio = dict.lookupOrDefault<scalar>(visualisationKeys::NORMALS_RATIO_KEY, defaultVisualisation->normalsRatio);
    normalsColor = dict.lookupOrDefault<point>(visualisationKeys::NORMALS_COLOR_KEY, defaultVisualisation->normalsColor);
    normalsOpacity = dict.lookupOrDefault<scalar>(visualisationKeys::NORMALS_OPACITY_KEY, defaultVisualisation->normalsOpacity);

    activeColour = dict.lookupOrDefault(visualisationKeys::ACTIVE_COLOR_KEY, defaultVisualisation->activeColour);
    showActiveGIBBoundary = dict.lookupOrDefault(visualisationKeys::SHOW_ACTIVE_KEY, defaultVisualisation->showActiveGIBBoundary);

    surfaceLic.readFromDict(dict, defaultVisualisation->surfaceLic);
}

const Visualisation* Visualisation::solidColorDefaults()
{
    static const Visualisation visualisation = createSolidColorDefaultValues();
    return &visualisation;
}

const Visualisation* Visualisation::indexDefaults()
{
    static const Visualisation visualisation = createIndexDefaultValues();
    return &visualisation;
}

bool Visualisation::isVisible() const
{
    return visible && opacity > 0.0;
}

bool Visualisation::isTransparent() const
{
    return isVisible() && (opacity < 1.0 || (showNormals && normalsOpacity < 1.0 && normalsOpacity > 0.0));
}

Visualisation& Visualisation::operator=(const Visualisation& other)
{
    if(this != &other)
    {
        this->visible = other.visible;
        this->showActiveGIBBoundary = other.showActiveGIBBoundary;
        this->opacity = other.opacity;
        this->colour = other.colour;
        this->activeColour = other.activeColour;
        this->lineColour = other.lineColour;
        this->lineThickness = other.lineThickness;
        this->colourField = other.colourField;
        this->supplementaryInfo = other.supplementaryInfo;
        this->backfaceStyling = other.backfaceStyling;
        this->surfaceLic = other.surfaceLic;
    }
    return *this;
}

Visualisation Visualisation::createCommonDefaultValues()
{
    Visualisation visualisation;
    visualisation.visible = false;
    visualisation.showActiveGIBBoundary = false;
    visualisation.opacity = 1.0;
    visualisation.colour = point(GRAY[0], GRAY[1], GRAY[2]);
    visualisation.activeColour = point(MAGENTA[0], MAGENTA[1], MAGENTA[2]);
    visualisation.lineColour = point(BLACK[0], BLACK[1], BLACK[2]);
    visualisation.lineThickness = 1.0;

    visualisation.showNormals = false;
    visualisation.normalsLength = 0.01;
    visualisation.normalsRatio = 1;
    visualisation.normalsOpacity = 1.0;
    visualisation.normalsColor = point(GRAY[0], GRAY[1], GRAY[2]);

    return {visualisation};
}

Visualisation Visualisation::createSolidColorDefaultValues()
{
    Visualisation visualisation = createCommonDefaultValues();
    visualisation.colourField = foamField(string("solidColor.SOLID"));
    return visualisation;
}

Visualisation Visualisation::createIndexDefaultValues()
{
    Visualisation visualisation = createCommonDefaultValues();
    visualisation.colourField = foamField(string("index.INDEXED"));
    return visualisation;
}

FoamFields Visualisation::getRequiredFields() const
{
    FoamFields fields;
    if (!colourField.isSolidColour())
    {
        fields.addField(colourField);
    }
    if (surfaceLic.visible)
    {
        fields.addField(foamField::withPointAssociation(surfaceLic.vectorField));
    }
    return fields;
}

} // End namespace Foam
