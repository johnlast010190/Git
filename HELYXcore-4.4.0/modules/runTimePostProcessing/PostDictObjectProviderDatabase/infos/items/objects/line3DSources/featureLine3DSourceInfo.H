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
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::FeatureLineSourceInfo

Description
    Contains the info for a 3D line from a feature line

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/objects/line3DSources/line3DSourceInfo.H"
#include "dataStructs/objects/line3D/featureLine3DSourceData.H"
#include "db/dictionary/dictionary.H"
#include "itemDataSetProviders/objects/line3d/featureLine3DDataSetProvider.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class FeatureLine3DSourceInfo : public Line3DSourceInfo
{
public:

    explicit FeatureLine3DSourceInfo(const dictionary& line3DDict) : data_(line3DDict) {};

    ~FeatureLine3DSourceInfo() override = default;

    void readDict(const dictionary& line3DDict, const string& caseFolder) override
    {
        data_.readDict(line3DDict);
    }

    void addExtraSources(std::vector<Id> &sources) const override
    {
        sources.push_back(data_.featureLineId);
    }

    void computeAndAddToHash(size_t& hash) const override
    {
        data_.computeAndAddToHash(hash);
    }

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider(const std::string& name) const override
    {
        return new FeatureLine3DDataSetProvider(name);
    }

    bool operator==(const Line3DSourceInfo& other) const override
    {
        const auto* casted = dynamic_cast<const FeatureLine3DSourceInfo*>(&other);
        if (casted == nullptr)
        {
            return false;
        }
        return data_ == casted->data_;
    }
private:
    FeatureLine3DSourceData data_;
};

} // End namespace Foam

// ************************************************************************* //
