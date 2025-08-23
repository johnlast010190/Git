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
    Foam::functionObjects::runTimeVis::StreamlineSourceInfo

Description
    Contains the streamline point source infos, based on a template class

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "itemDataSetProviders/objects/streamlineSources/streamlineSourceProviders.H"
#include "db/dictionary/dictionary.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class TemplatedStreamlineSourceInfo
{
public:
    virtual ~TemplatedStreamlineSourceInfo() = default;


    // Member Functions

    [[nodiscard]] virtual std::shared_ptr<StreamlineSourceProvider> createStreamlineSourceProvider() const = 0;

    virtual void computeAndAddToHash(size_t& hash) const  = 0;

    virtual bool isDataEqualTo(const TemplatedStreamlineSourceInfo* other) const = 0;

    virtual void addToSources(std::vector<Id>& vtkNotUsed(sources)) const {};
};

template <typename Data, typename Provider>
class TemplateStreamlineSourceInfo : public TemplatedStreamlineSourceInfo
{
protected:
    Data data_;

public:
    explicit TemplateStreamlineSourceInfo(const dictionary& streamlinesDict)
    :
        data_(streamlinesDict)
    {};

    [[nodiscard]] std::shared_ptr<StreamlineSourceProvider> createStreamlineSourceProvider() const override
    {
        return std::shared_ptr<StreamlineSourceProvider>(new Provider(data_));
    };

    void computeAndAddToHash(size_t& hash) const override
    {
        data_.computeAndAddToHash(hash);
    }

    bool isDataEqualTo(const TemplatedStreamlineSourceInfo* other) const override
    {
        auto casted = dynamic_cast<const TemplateStreamlineSourceInfo<Data, Provider>*>(other);
        if (casted == nullptr)
        {
            return false;
        }
        return data_ == casted->data_;
    };
};

// Declaration of the individual cutting infos
typedef TemplateStreamlineSourceInfo<LineSourceData, LineSourceProvider> LineSourceInfo;
typedef TemplateStreamlineSourceInfo<PointCloudSourceData, PointCloudSourceProvider> PointCloudSourceInfo;

class PatchSourceInfo : public TemplateStreamlineSourceInfo<PatchSourceData, PatchSourceProvider>
{
public:
    explicit PatchSourceInfo(const dictionary& streamlinesDict)
    :
        TemplateStreamlineSourceInfo<PatchSourceData, PatchSourceProvider>(streamlinesDict)
    {};

    void addToSources(std::vector<Id>& sources) const override
    {
        sources.push_back(data_.patchId);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
