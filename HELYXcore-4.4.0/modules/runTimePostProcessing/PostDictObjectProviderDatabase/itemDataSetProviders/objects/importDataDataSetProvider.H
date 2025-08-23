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
    Foam::functionObjects::runTimeVis::ImportDataDataSetProvider

Description
    Provider for the threshold 3d objects, which computes their vtk dataset

SourceFiles
    importDataDataSetProvider.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "itemDataSetProviders/externalItemDataSetProvider.H"

#include "postDict/postDictKeys.H"
#include "db/dictionary/dictionary.H"
#include "hash/hasher.H"

#include "engysExternalDataSetReader.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct ImportDataData
{
    engysExternalDataSetFileFormat format;
    fileName path;
    scalar timeStep;

    void readDict(const dictionary& objectDict)
    {
        word formatKey = objectDict.lookupOrDefault<word>(importDataObjectKeys::FORMAT_KEY, importDataObjectKeys::PVD_KEY);

        if (formatKey == importDataObjectKeys::PVD_KEY || formatKey == importDataObjectKeys::FIELD_SAMPLING_KEY)
        {
            format = engysExternalDataSetFileFormat::PVD;
        }
        else if (formatKey == importDataObjectKeys::VTK_KEY)
        {
            format = engysExternalDataSetFileFormat::VTK;
        }
        else if (formatKey == importDataObjectKeys::ENSIGHT_KEY)
        {
            format = engysExternalDataSetFileFormat::ENSIGHT; 
        }
        
        path = objectDict.lookup<fileName>(importDataObjectKeys::PATH_KEY);
        timeStep = objectDict.lookupOrDefault<scalar>(importDataObjectKeys::TIME_STEP_KEY, -1.0);
    }

    bool operator==(const ImportDataData& other) const
    {
        return format == other.format && path == other.path && timeStep == other.timeStep;
    }

    [[nodiscard]] size_t computeAndReturnHash() const
    {
        size_t hash = 0;
        hasher::hash_combine(hash, format);
        hasher::hash_combine(hash, path);
        hasher::hash_combine(hash, timeStep);
        return hash;
    }
};

/*---------------------------------------------------------------------------*\
                 Class ImportDataDataSetProvider Declaration
\*---------------------------------------------------------------------------*/

class ImportDataDataSetProvider : public ExternalItemDataSetProvider
{
private:
    vtkSmartPointer<engysExternalDataSetReader> reader_;

protected:
        void update(scalar currentTime) override;

public:
        ImportDataDataSetProvider(const std::string& name, const ImportDataData& dictData);
        ~ImportDataDataSetProvider() override = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
