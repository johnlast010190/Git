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
    (c) 2019 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

File
    Foam::functionObjects::runTimeVis::hasher

Description
    Contains a function used to calculate and combine hashes

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include <string>

#include "primitives/Vector/vector/vector.H"

namespace Foam::functionObjects::runTimeVis::hasher
{

    // The constant 0x9e3779b9 comes from boost
#define CONSTANT 0x9e3779b9

    inline void hash_combine(std::size_t& seed, const std::string& v) {
        std::hash<std::string> hasher;
        seed ^= hasher(v) + CONSTANT + (seed<<6) + (seed>>2);
    }

    inline void hash_combine(std::size_t& seed, double v) {
        std::hash<double> hasher;
        seed ^= hasher(v) + CONSTANT + (seed<<6) + (seed>>2);
    }

    inline void hash_combine(std::size_t& seed, float v) {
        std::hash<float> hasher;
        seed ^= hasher(v) + CONSTANT + (seed<<6) + (seed>>2);
    }

    inline void hash_combine(std::size_t& seed, long v) {
        std::hash<long> hasher;
        seed ^= hasher(v) + CONSTANT + (seed<<6) + (seed>>2);
    }

    inline void hash_combine(std::size_t& seed, int v) {
        std::hash<int> hasher;
        seed ^= hasher(v) + CONSTANT + (seed<<6) + (seed>>2);
    }

    inline void hash_combine(std::size_t& seed, vector v) {
        hash_combine(seed, v.x());
        hash_combine(seed, v.y());
        hash_combine(seed, v.z());
    }

    inline void hash_combine(std::size_t& seed, bool v) {
        std::hash<bool> hasher;
        seed ^= hasher(v) + CONSTANT + (seed<<6) + (seed>>2);
    }


    // These combine with a simple add, because then the order doesn't matter

    inline void hash_add_unordered(std::size_t& seed, double v) {
        std::hash<double> hasher;
        seed += hasher(v);
    }

    inline void hash_add_unordered(std::size_t& seed, std::size_t v) {
        seed += v;
    }

    inline void hash_add_unordered(std::size_t& seed, const std::string& v) {
        std::hash<std::string> hasher;
        seed += hasher(v);
    }
}
// End namespace
