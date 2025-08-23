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
    (c) 2022-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "colourLookupTablesInfo.H"
#include "postDict/postDictKeys.H"
#include "dictionaries/dictionaries.H"

#undef Log
#include "vtkLogger.h"

namespace Foam::functionObjects::runTimeVis
{

ColourLookupTablesInfo::ColourLookupTablesInfo
(
        const dictionary& dict
) : defaultLut_()
{
    if (dict.found(sceneKeys::COLOR_LEGEND_OPTIONS_DICT_KEY))
    {
        const dictionary& defaultLegendsDict = dict.subDict(sceneKeys::COLOR_LEGEND_OPTIONS_DICT_KEY);
        foamField defaultName(std::string("default"));
        defaultLut_ = ColourLookupTableData(defaultName, defaultLegendsDict);
    }

    const dictionary& legendsDict = GET_OPTIONAL_DICTIONARY(dict, sceneKeys::COLOR_LEGENDS_DICT_KEY);

    colourLutData_.clear();
    for(const std::string& fieldName : legendsDict.toc())
    {
        foamField foamFieldName(fieldName);
        ColourLookupTableData tempLut
        (
            foamField(foamFieldName.lessAssociation()),
            legendsDict.subDict(fieldName),
            defaultLut_
        );
        colourLutData_.insert(foamFieldName.lessAssociation(), tempLut);
    }
}

ColourLookupTableData ColourLookupTablesInfo::getColorLookupTableData
(
    const foamField& field,
    const ColourLookupTablesInfo& baseFields
) const
{
    word fieldLA = field.lessAssociation();
    if (colourLutData_.found(fieldLA))
    {
        return colourLutData_[fieldLA];
    }
    else
    {
        return baseFields.getColorLookupTableData(fieldLA);
    }
}

ColourLookupTableData ColourLookupTablesInfo::getColorLookupTableData(const word& fieldLA) const
{
    if (colourLutData_.found(fieldLA))
    {
        return colourLutData_[fieldLA];
    }
    else
    {
        vtkLogF(INFO, "Creating from default lut for %s", fieldLA.c_str());
        foamField foamFieldName(fieldLA);
        ColourLookupTableData tempLut(foamFieldName, {}, defaultLut_);
        return tempLut;
    }
}

}
