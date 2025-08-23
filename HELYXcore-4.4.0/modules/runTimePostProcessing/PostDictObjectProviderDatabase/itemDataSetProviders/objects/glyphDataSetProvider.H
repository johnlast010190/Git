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
    (c) 2022 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::GlyphDataSetProvider

Description
    Provider for the vectors 3d objects, which computes their vtk dataset

SourceFiles
    glyphDataSetProvider.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "itemDataSetProviders/itemDataSetProvider.H"
#include "types/glyphType.H"
#include "types/glyphSourceType.H"

#include "hash/hasher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class engysGlyphs;

namespace Foam::functionObjects::runTimeVis
{

struct GlyphObjectData
{
    GlyphType glyphType;
    foamField orientationField;
    bool enableScalingByField;
    foamField scalingField;
    scalar maxGlyphLength;
    struct GlyphDistribution
    {
        GlyphSourceType type;
        label maximumNumberOfPoints = 0;
        label seed = 0;
        label pointsInterval = 0;

        void readDict(const dictionary& glyphDistDict)
        {
            type = glyphDistDict.lookup<GlyphSourceType>("type");
            switch(type.getValue())
            {
            case GlyphSourceType::UNIFORM_BOUNDS_BASED:
                maximumNumberOfPoints = glyphDistDict.lookup<label>("maximumNumberOfPoints");
                seed = glyphDistDict.lookup<label>("seed");
                pointsInterval = 0;
                break;
            case GlyphSourceType::EVERY_NTH_POINT:
                maximumNumberOfPoints = 0;
                seed = 0;
                pointsInterval = glyphDistDict.lookup<label>("pointsInterval");
                break;
            case GlyphSourceType::ALL_POINTS:
                maximumNumberOfPoints = 0;
                seed = 0;
                pointsInterval = 0;
                break;
            }
        }

        bool operator==(const GlyphDistribution& other) const
        {
            if (type != other.type)
            {
                return false;
            }
            switch(type.getValue())
            {
            case GlyphSourceType::UNIFORM_BOUNDS_BASED:
                return maximumNumberOfPoints == other.maximumNumberOfPoints &&
                                        seed == other.seed;
            case GlyphSourceType::EVERY_NTH_POINT:
                return pointsInterval == other.pointsInterval;
            case GlyphSourceType::ALL_POINTS:
            default:
                return true;
            }
        }
        void computeAndAddToHash(size_t& hash) const
        {
            hasher::hash_combine(hash, type.getValue());
            switch(type.getValue())
            {
            case GlyphSourceType::UNIFORM_BOUNDS_BASED:
                hasher::hash_combine(hash, maximumNumberOfPoints);
                hasher::hash_combine(hash, seed);
                break;
            case GlyphSourceType::EVERY_NTH_POINT:
                hasher::hash_combine(hash, pointsInterval);
                break;
            case GlyphSourceType::ALL_POINTS:
                break;
            }
        }
    } glyphDistribution;

    void readDict(const dictionary& dict)
    {
        glyphType = dict.lookup<GlyphType>("glyphType");
        orientationField = foamField(dict.lookup<word>("orientationField"));
        enableScalingByField = dict.lookup<bool>("enableScalingByField");
        scalingField = foamField(dict.lookup<string>("scalingField"));
        maxGlyphLength = dict.lookup<scalar>("maxGlyphLength");

        const dictionary& glyphDistDict = dict.subDict("glyphDistribution");
        glyphDistribution.readDict(glyphDistDict);
    }

    bool operator==(const GlyphObjectData& other) const
    {
        return          glyphType == other.glyphType &&
                 orientationField == other.orientationField &&
             enableScalingByField == other.enableScalingByField &&
                     scalingField == other.scalingField &&
                   maxGlyphLength == other.maxGlyphLength &&
                glyphDistribution == other.glyphDistribution;

    }
    [[nodiscard]] size_t computeAndReturnHash() const
    {
        size_t hash = 0;
        hasher::hash_combine(hash, glyphType.getValue());
        hasher::hash_combine(hash, orientationField);
        hasher::hash_combine(hash, enableScalingByField);
        hasher::hash_combine(hash, scalingField);
        hasher::hash_combine(hash, maxGlyphLength);
        glyphDistribution.computeAndAddToHash(hash);

        return hash;
    }
};

/*---------------------------------------------------------------------------*\
                     Class GlyphDataSetProvider Declaration
\*---------------------------------------------------------------------------*/

class GlyphDataSetProvider : public ItemDataSetProvider
{
private:

protected:

        vtkSmartPointer<engysGlyphs> vectorsFilter_;

        GlyphObjectData glyphData_;

        void update(scalar currentTime) override;

public:
        GlyphDataSetProvider(const std::string& name, const GlyphObjectData& dictData);
        ~GlyphDataSetProvider() override = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimeVis
// End namespace functionObjects
// End namespace Foam

// ************************************************************************* //
