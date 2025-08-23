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
    Foam::functionObjects::runTimeVis::FeatureLineData

Description
    Contains the dict data necessary for feature lines

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "hash/hasher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations

namespace Foam::functionObjects::runTimeVis
{

struct FeatureLineData {
    string nameWithExtension;
    string caseFolder;

    FeatureLineData(const string& name, const string& caseFolder)
    {
        this->nameWithExtension = name + ".eMesh";
        this->caseFolder = caseFolder;
    }

    bool operator==(const FeatureLineData& other) const
    {
        return nameWithExtension == other.nameWithExtension &&
                      caseFolder == other.caseFolder;
    }

    [[nodiscard]] size_t computeAndReturnHash() const
    {
        size_t hash = 0;
        hasher::hash_combine(hash, nameWithExtension);
        hasher::hash_combine(hash, caseFolder);
        return hash;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
