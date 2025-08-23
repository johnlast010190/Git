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

Class
    Foam::functionObjects::runTimeVis::FieldSamplingObjectData

Description
    Contains the dict data for the field sampling

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "postDict/postDictKeys.H"
#include "db/dictionary/dictionary.H"
#include "hash/hasher.H"
#include "types/fieldSamplingToleranceType.H"

#include "meshes/boundBox/boundBox.H"
#include "meshes/primitiveShapes/point/point.H"
#include "primitives/Vector/vector/vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct FieldSamplingObjectData
{
    string caseFolder;
    boundBox samplingBounds;
    List<word> desiredFields;
    label elementsX = 0;
    label elementsY = 0;
    label elementsZ = 0;
    FieldSamplingToleranceType toleranceType;
    double toleranceValue;

    void readDict(const dictionary &fieldSamplingDict)
    {
        samplingBounds.clear();
        samplingBounds.add(fieldSamplingDict.lookup<point>(fieldSamplingObjectKeys::SAMPLING_BOUNDS_MIN_KEY));
        samplingBounds.add(fieldSamplingDict.lookup<point>(fieldSamplingObjectKeys::SAMPLING_BOUNDS_MAX_KEY));
        auto elements = fieldSamplingDict.lookup<Vector<label>>(fieldSamplingObjectKeys::ELEMENTS_KEY);
        elementsX = elements.x() + 1;
        elementsY = elements.y() + 1;
        elementsZ = elements.z() + 1;
        if (fieldSamplingDict.found(fieldSamplingObjectKeys::FIELDS_KEY))
        {
            desiredFields = fieldSamplingDict.lookup<List<word>>(fieldSamplingObjectKeys::FIELDS_KEY);
        }
        else
        {
            desiredFields.clear();
        }
        if (fieldSamplingDict.found(fieldSamplingObjectKeys::TOLERANCE_TYPE))
        {
            toleranceType = fieldSamplingDict.lookup<FieldSamplingToleranceType>(fieldSamplingObjectKeys::TOLERANCE_TYPE);
        }
        else
        {
            bool useTolerance = fieldSamplingDict.lookupOrDefault(fieldSamplingObjectKeys::USE_TOLERANCE, false);
            if (useTolerance) toleranceType = FieldSamplingToleranceType(FieldSamplingToleranceType::Value::ABSOLUTE);
            else toleranceType = FieldSamplingToleranceType(FieldSamplingToleranceType::Value::OFF);
        }
        toleranceValue = fieldSamplingDict.lookupOrDefault(fieldSamplingObjectKeys::TOLERANCE_VALUE, 0.0);
    }

    bool operator==(const FieldSamplingObjectData &other) const
    {
        return samplingBounds == other.samplingBounds &&
               elementsX == other.elementsX &&
               elementsY == other.elementsY &&
               elementsZ == other.elementsZ &&
               toleranceType == other.toleranceType &&
               toleranceValue == other.toleranceValue;
    }

    void computeAndAddToHash(size_t &hash) const
    {
        hasher::hash_combine(hash, samplingBounds.min());
        hasher::hash_combine(hash, samplingBounds.max());
        hasher::hash_combine(hash, elementsX);
        hasher::hash_combine(hash, elementsY);
        hasher::hash_combine(hash, elementsZ);
        hasher::hash_combine(hash, toleranceType.getValue());
        hasher::hash_combine(hash, toleranceValue);
    }
};

} // End namespace

// ************************************************************************* //
