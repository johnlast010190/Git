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
    (c) 2023-2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::TurboSliceSpanwiseObjectData

Description
    Contains the dict data for the spanwise slice plot post-processing object

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "postDict/postDictKeys.H"
#include "db/dictionary/dictionary.H"
#include "hash/hasher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct TurboSliceSpanwiseObjectData
{
    scalar streamSlice;
    bool unravel;

    void readDict(const dictionary &fieldSamplingDict)
    {
        streamSlice = fieldSamplingDict.lookup<scalar>(turbopropKeys::STREAM_VALUE_KEY);
        unravel = fieldSamplingDict.lookup<bool>(turbopropKeys::BLADE_TO_BLADE_KEY);
    }

    bool operator==(const TurboSliceSpanwiseObjectData &other) const
    {
        return streamSlice == other.streamSlice &&
               unravel == other.unravel;
    }

    void computeAndAddToHash(size_t &hash) const
    {
        hasher::hash_combine(hash, streamSlice);
        hasher::hash_combine(hash, unravel);
    }
};

} // End namespace Foam

// ************************************************************************* //
