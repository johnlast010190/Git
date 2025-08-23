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
    Foam::functionObjects::runTimeVis::FontInfo

Description
    Contains the font configuration data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimeVis_fontInfo_H
#define functionObjects_runTimeVis_fontInfo_H

#include "postDict/postDictKeys.H"
#include "types/font.H"
#include "primitives/Vector/Vector.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct FontInfo
{
    Fonts type;
    label size = 12;
    scalar opacity = 1.0;
    Vector<scalar> colour = {0,0,0};
    bool bold = false;
    bool italics = false;
    bool shadow = false;

    FontInfo() = default;

    void readDict
    (
        const dictionary& fontDict,
        const FontInfo& defaultFont
    )
    {
        type = fontDict.lookupOrDefault(fontKeys::FONT_KEY, defaultFont.type);
        size = fontDict.lookupOrDefault(fontKeys::SIZE_KEY, defaultFont.size);
        colour = fontDict.lookupOrDefault(fontKeys::COLOR_KEY, defaultFont.colour);
        opacity = fontDict.lookupOrDefault(fontKeys::OPACITY_KEY, defaultFont.opacity);
        bold = fontDict.lookupOrDefault(fontKeys::BOLD_KEY, defaultFont.bold);
        italics = fontDict.lookupOrDefault(fontKeys::ITALICS_KEY, defaultFont.italics);
        shadow = fontDict.lookupOrDefault(fontKeys::SHADOW_KEY, defaultFont.shadow);
    }
};

} // End namespace

#endif
