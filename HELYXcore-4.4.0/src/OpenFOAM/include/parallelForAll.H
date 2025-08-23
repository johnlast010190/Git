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
    (c) 2021 Engys Ltd.

Namespace
    parallelForAll

Description
    Tool for code parallelization with TBB

\*---------------------------------------------------------------------------*/

#ifndef ForEach_H
#define ForEach_H

#include "stdFoam.H"

#ifdef HELYX_USE_TBB
    #include <tbb/parallel_for.h>
    #include <tbb/parallel_for_each.h>
    #include <tbb/blocked_range.h>
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    template<typename T, typename Func>
    void parallelForAll
    (
        const UList<T>& list,
        const Func& func
    )
    {
#ifndef HELYX_USE_TBB
        forAll(list, i)
        {
            func(i);
        }
#else
        try
        {
            tbb::parallel_for
            (
                label(0),
                list.size(),
                func
            );
        }
        catch (std::exception& ex)
        {
            throw error(ex.what());
        }
#endif
    }
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
