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
    (c) 2022-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "renderInfo.H"
#include "postDict/postDictKeys.H"
#include "types/exportType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

RenderInfo::RenderInfo(const dictionary& postDict)
{
    auto exportFormats = postDict.lookupOrDefault<List<ExportType>>(
        controlDictKeys::EXPORT_FORMATS_KEY,
        List<ExportType>{ExportType(ExportType::PNG)}
    );
    for(const ExportType& exportFormat : exportFormats)
    {
        switch (exportFormat.getValue())
        {
            case ExportType::PNG:
                exportPng = true;
                // Read image width and height
                width = postDict.lookup<label>(imageExportKeys::IMAGE_WIDTH_KEY);
                height = postDict.lookup<label>(imageExportKeys::IMAGE_HEIGHT_KEY);

                transparentBackground = postDict.lookupOrDefault(imageExportKeys::TRANSPARENT_KEY, false);
                cropToContents = postDict.lookupOrDefault(imageExportKeys::CROP_TO_CONTENTS_KEY, false);
                centerToCroppedContents = postDict.lookupOrDefault(
                    imageExportKeys::CENTER_CROPPED_CONTENTS_KEY,
                    false
                );
                break;
            case ExportType::X3D:
                exportX3d = true;
                break;
            case ExportType::EDF:
                exportEdf = true;
                width = postDict.lookup<label>(imageExportKeys::IMAGE_WIDTH_KEY);
                height = postDict.lookup<label>(imageExportKeys::IMAGE_HEIGHT_KEY);
                if (postDict.found(imageExportKeys::EDF_FIELDS_KEY))
                {
                    edfFields = postDict.lookup<List<foamField>>(imageExportKeys::EDF_FIELDS_KEY);
                }
                edfCompressionLevel = EDFCompressionLevel(
                    postDict.lookupOrDefault<int>(
                        imageExportKeys::EDF_COMPRESSION_LEVEL_KEY,
                        BALANCED
                    ));
                break;
            case ExportType::UNKNOWN:
                break;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
