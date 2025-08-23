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

#include "UPstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(UPstream, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::UPstream::commsTypes,
        3
    >::names[] =
    {
        "blocking",
        "scheduled",
        "nonBlocking"
    };
}

const Foam::NamedEnum<Foam::UPstream::commsTypes, 3>
    Foam::UPstream::commsTypeNames;


bool Foam::UPstream::parRun_(false);

bool Foam::UPstream::haveThreads_(false);

Foam::LIFOStack<Foam::label> Foam::UPstream::freeComms_;

Foam::DynamicList<int> Foam::UPstream::myProcNo_(10);

Foam::DynamicList<Foam::List<int>> Foam::UPstream::procIDs_(10);

Foam::DynamicList<Foam::label> Foam::UPstream::parentCommunicator_(10);

int Foam::UPstream::msgType_(1);


Foam::DynamicList<Foam::List<Foam::UPstream::commsStruct>>
Foam::UPstream::linearCommunication_(10);

Foam::DynamicList<Foam::List<Foam::UPstream::commsStruct>>
Foam::UPstream::treeCommunication_(10);


// Allocate a serial communicator. This gets overwritten in parallel mode
// (by UPstream::setParRun())
Foam::UPstream::communicator serialComm
(
    -1,
    Foam::labelList(1, Foam::Zero),
    false
);


bool Foam::UPstream::floatTransfer
(
    Foam::debug::optimisationSwitch("floatTransfer", 0)
);
registerOptSwitch
(
    "floatTransfer",
    bool,
    Foam::UPstream::floatTransfer
);

int Foam::UPstream::nProcsSimpleSum
(
    Foam::debug::optimisationSwitch("nProcsSimpleSum", 16)
);
registerOptSwitch
(
    "nProcsSimpleSum",
    int,
    Foam::UPstream::nProcsSimpleSum
);

Foam::UPstream::commsTypes Foam::UPstream::defaultCommsType
(
    commsTypeNames.read(Foam::debug::optimisationSwitches().lookup("commsType"))
);
