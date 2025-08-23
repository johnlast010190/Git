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

\*---------------------------------------------------------------------------*/

#ifndef TBBTimer_H
#define TBBTimer_H

#include "cpuTime/cpuTime.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include <tbb/tick_count.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    struct Timer
    {
        tbb::tick_count start_, end_;
        scalar duration_;
        cpuTime cpuTimer_; //https://pubs.opengroup.org/onlinepubs/9699919799/functions/times.html
        word text_;

        Timer(word text="")
        {
            start_ = tbb::tick_count::now();
            cpuTimer_ = cpuTime();
            text_ = text;
        }

        ~Timer()
        {
            end_ = tbb::tick_count::now();
            duration_ = (end_ - start_).seconds();

            scalar cpuTime = cpuTimer_.elapsedCpuTime();

            Info<< "=== Timer: "
                << setw(38) << text_
                << ", elapsed " << setw(9) << duration_ << setprecision(3) << " s "
                << "(cpu " << setw(5) << cpuTime << " s)"
                << endl;
        }
    };
}
 // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
