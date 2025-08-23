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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#ifndef HELYXCORE_MPIRequest_H
#define HELYXCORE_MPIRequest_H

#include <mpi.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/**
 * Represents a result value that will eventually exist for an MPI request.
 * @tparam T
 */
class MPIRequest {
    MPI_Request req{MPI_REQUEST_NULL};

public:

    /// Wait for the underlying request to finish, and return the result value.
    /// This may not be called from multiple threads, although that could be changed.
    /// This may not be called twice.
    void await() {
        if (req == MPI_REQUEST_NULL) {
            return;
        } else {
            MPI_Wait(&req, MPI_STATUS_IGNORE);
        }
    }

    /**
     * Checks if the request has finished yet or not. Allows you to test if the value
     * is opportunistically ready to read without having to block.
     *
     * Use of this function is discouraged, since it is inherently a race condition. It is
     * fairly easy to accidentally introduce nondeterminism into your program if you use
     * this function. It is also possible to make things faster, if done correctly :D.
     */
    bool isReady() {
        int out;
        MPI_Request_get_status(req, &out, MPI_STATUS_IGNORE);
        return out;
    }

    operator MPI_Request*() {
        return &req;
    }
};



#endif //HELYXCORE_MPIRequest_H
