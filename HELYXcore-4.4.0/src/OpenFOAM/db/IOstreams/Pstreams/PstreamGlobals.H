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
    (c) 2013-2015 OpenFOAM Foundation

Namespace
    Foam::PstreamGlobals

Description
    Global functions and variables for working with parallel streams,
    but principally for mpi

SourceFiles
    PstreamGlobals.C

\*---------------------------------------------------------------------------*/

#ifndef PstreamGlobals_H
#define PstreamGlobals_H

#include "containers/Lists/DynamicList/DynamicList.H"

#include <mpi.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PstreamGlobals Declaration
\*---------------------------------------------------------------------------*/

namespace PstreamGlobals
{

extern DynamicList<MPI_Request> outstandingRequests_;

//extern int nRequests_;
//extern DynamicList<label> freedRequests_;

extern int nTags_;

extern DynamicList<int> freedTags_;


// Current communicators. First element will be MPI_COMM_WORLD
extern DynamicList<MPI_Comm> MPICommunicators_;
extern DynamicList<MPI_Group> MPIGroups_;

void checkCommunicator(const label, const label procNo);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
