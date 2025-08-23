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
    (c) 2015-2019 OpenCFD Ltd.
    (c) 2024-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "pvdInfo.H"

#include "postDict/postDictKeys.H"
#include "types/exportType.H"
#include "db/dictionary/dictionary.H"
#include "infos/sceneInfo.H"

namespace Foam::functionObjects::runTimeVis
{

PvdInfo::PvdInfo(const dictionary& postDict)
{
    auto exportFormats = postDict.lookupOrDefault<List<ExportType>>(
        controlDictKeys::EXPORT_FORMATS_KEY,
        List<ExportType>{ExportType(ExportType::UNKNOWN)}
    );
    exportPvd = exportFormats.found(ExportType(ExportType::PVD));
    if (exportPvd)
    {
        pvdFieldsType = postDict.lookupOrDefault(pvdExportKeys::PVD_FIELDS_TYPE_KEY, PvdFieldsType());
        writePointData = postDict.lookupOrDefault(pvdExportKeys::PVD_WRITE_POINT_DATA_KEY, true);
        oneFilePerProcess = postDict.lookupOrDefault(pvdExportKeys::PVD_WRITE_ONE_FILE_PER_PROCESS_KEY, false);
        removeOldObsoleteFiles = postDict.lookupOrDefault(pvdExportKeys::PVD_REMOVE_OLD_OBSOLETE_FILES_KEY, false);
        switch (pvdFieldsType.getValue())
        {
            case PvdFieldsType::ALL:
            {
                pvdRequirements.addToRequiredFields(foamField::AllFieldsMarker());
                break;
            }
            case PvdFieldsType::FROM_LIST:
            {
                if (!postDict.found(pvdExportKeys::PVD_FIELDS_LIST_KEY))
                {
                    FatalError << "PVD set to export fields specified from a list, but no list was defined "
                               << "(" << pvdExportKeys::PVD_FIELDS_LIST_KEY << ")"
                               << endl << abort(FatalError);
                }
                auto pvdFieldsList = postDict.lookup<List<word>>(pvdExportKeys::PVD_FIELDS_LIST_KEY);
                for (const word &pvdField: pvdFieldsList)
                {
                    pvdRequirements.addToRequiredFields(foamField(pvdField));
                }
                break;
            }
            case PvdFieldsType::FROM_SCENE:
                // No need to add any extra required fields
                break;
        }
        compressionLevel = postDict.lookupOrDefault(
            pvdExportKeys::PVD_COMPRESSION_LEVEL_KEY,
            CompressionLevelType(CompressionLevelType::BALANCED)
            );
    }
}

const ItemRequirements& PvdInfo::getItemRequirements() const
{
    return pvdRequirements;
}

FoamFields PvdInfo::getAllowedFieldsForPVD(const SceneInfo& sceneInfo) const
{
    if (exportPvd)
    {
        switch (pvdFieldsType.getValue())
        {
            case PvdFieldsType::ALL:
                return allFields();
            case PvdFieldsType::FROM_LIST:
                return this->pvdRequirements.getRequiredFields();
            case PvdFieldsType::FROM_SCENE:
               return sceneInfo.getItemRequirements().getRequiredFields();
        }
    }
    return {};
}

FoamFields createWithAllFields()
{
    FoamFields allFields;
    allFields.addAllFieldsMarker();
    return allFields;
}

const FoamFields& PvdInfo::allFields()
{
    static FoamFields allFields = createWithAllFields();
    return allFields;
}

} // End namespace Foam