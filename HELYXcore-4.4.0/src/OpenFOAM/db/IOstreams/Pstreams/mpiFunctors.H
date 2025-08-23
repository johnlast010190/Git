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

#ifndef MPI_FUNCTORS_H
#define MPI_FUNCTORS_H

// Trampoline functions (as functors) for calling various MPI functions.
//
// Makes it possible to pass equivalently-shaped MPI primitives as template parmaeters, which can reduce
// duplication a lot. Also serves as a layer of indirection between this library and MPI to babysit the insanity
// of the dynamic-code mechanism (which does not have reasonable include management).
//
// These would be better if defined as templates, but the "dynamic code" system's inability to handle mpi.h prevents this.
// Were that done, we could get some improved type-safety and TBAA properties out of this, and a teeny-tiny speedup:
// this version creates at least one redundant stack frame per call.


namespace Foam {

struct mpiAllGatherv {
    int operator()(
        const void* sendbuf, int sendcount,
        void* recvbuf, const int recvcounts[],
        const int displs[],
        int comm
    );
};

struct mpiAllGather {
    int operator()(
        const void* sendbuf, int  sendcount,
        void* recvbuf, int recvcount,
        int comm
    );
};

struct mpiGatherv {
    int operator()(
        const void* sendbuf, int sendcount,
        void* recvbuf, const int recvcounts[],
        const int displs[],
        int comm
    );
};

struct mpiGather {
    int operator()(
        const void* sendbuf, int sendcount,
        void* recvbuf, int recvcount,
        int comm
    );
};

} // Namespace Foam

#endif // MPI_FUNCTORS_H
