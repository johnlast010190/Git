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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.

Description
    Exchange data.

\*---------------------------------------------------------------------------*/

#include "db/IOstreams/Pstreams/PstreamReduceOps.H"
#include <iostream>

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace {

template<typename T>
constexpr T divRoundUp(T x, T y) {
    return (x + y - 1) / y;
}

template<int HowMany>
struct SomeBytes {
    std::byte bytes[HowMany];
};

} // Anonymous namespace

template<class Container, class T>
void Foam::Pstream::exchange
(
    const UList<Container>& sendBufs,
    const labelUList& recvSizes,
    List<Container>& recvBufs,
    label comm
)
{
    static_assert(contiguous<T>(), "exchange() only supports memcpyable payload types.");
    static_assert(std::is_base_of_v<Foam::UList<T>, Container>, "exchange() only supports list-esque containers.");

    int numSendBufs = sendBufs.size();
    if (numSendBufs != UPstream::nProcs(comm))
    {
        FatalErrorInFunction
                << "Size of list " << numSendBufs
                << " does not equal the number of processors "
                << UPstream::nProcs(comm)
                << Foam::abort(FatalError);
    }
    recvBufs.setSize(numSendBufs);

    if (!UPstream::parRun() || UPstream::nProcs(comm) == 1)
    {
        // If there's only one process, there is very little to do.
        recvBufs[Pstream::myProcNo(comm)] = sendBufs[Pstream::myProcNo(comm)];
        return;
    }

    if (UPstream::useMpiAllToAll)
    {
        exchangeAllToAll<Container, T>(sendBufs, recvSizes, recvBufs, comm);
    }
    else
    {
        exchangePointToPoint<Container, T>(sendBufs, recvSizes, recvBufs, comm);
    }
}

template<class Container, class T>
void Foam::Pstream::exchangeAllToAll
(
    const UList<Container>& sendBufs,
    const labelUList& recvSizes,
    List<Container>& recvBufs,
    label comm
)
{
    static_assert(sizeof(T) <= MaxMpiCustomTypeLengthInBytes, "That type is too long for an MPI exchange. Increase MaxMpiCustomTypeLengthInBytes.");

    int numSendBufs = sendBufs.size();

    // Since the payload type is contiguous, and the underlying container has a memcpyable
    // data buffer, all we actually need to do is shovel the bytes from A to B, which MPI
    // "conveniently" provides a function for.
    //
    // If you need to do stuff like exchange<HashMap<T>>, it's _probably_ most efficient
    // to list-ify it, exchange, and then rebuild the map on the receiver: the wasteful
    // reallocations etc. are likely cheaper than the multiple rounds of communication
    // that'd result from doing it incrementally, unless only a small fraction of a large
    // map is actually being transmitted.
    // `MPI_Alltoallw` may allow a smarter solution to that...
    //
    // Unfortunately, the MPI guys decided not to allow i64 length parameters, so to send
    // large messages one must use a datatype larger than 1 byte, and then send up to 2^32
    // of them (implying padding. Hooray).
    //
    // Below, "record" refers to "block of bytes used to cheat the length limit", and "element"
    // to "one of the things you're actually trying to send". We arrange for records to be
    // an integer multiple of element length for sanity. Increasing record length increases the
    // maximum possible message size as well as the maximum possible amount of wasted data
    // transmission (currently 127 bytes per processor, which is likely to be negligible, and
    // certainly less significant than the benefits of being able to leverage IP multicast and
    // friends to perform the operation). For a record length of `N` bytes, the maximum possible
    // message length is `2^(32 + N)` bytes.
    //
    // At present, almost all calls to this function come via PStreamBuffers, which just uses
    // it to schlep chars around, hiding the underlying datatype (which is almost always larger
    // than 1 byte). Were the true type being transmitted to be known to this function, you'd
    // effectively divide the necessary amount of wasted transmission by 2^(sizeof(T)), since you
    // could just inform MPI about T, and send records composed thereof.
    //
    // It may be possible to eliminate all of this via a binding to the FORTRAN API, but
    // 64-bit length inputs in the FORTRAN interfaces of MPI functions seem to be inconsistently
    // supported between implementations, and I was not issued a hazmat suit.
    //
    // </rant>

    // The maximum number of repetitions of T that fits into the maximum defined record type.
    constexpr label EltsPerRecord = MaxMpiCustomTypeLengthInBytes / sizeof(T);
    constexpr uint64_t RecordSize = sizeof(T) * EltsPerRecord;
    auto MpiRecordTy = UPstream::dataTypes[RecordSize];
    using RecordTy = SomeBytes<RecordSize>;

    label numRecordsToSend = 0;
    label numRecordsToRecv = 0;

    // TODO: Since each communicator can only do a single collective at a time, the buffer for all
    //       these counters could just be retained permanently, once-per-communicator, to avoid doing
    //       allocations every time we want to do a communication op.
    // For now just be moderately clever and do only one heap allocation for the counters.
    Foam::List<int> variousIntegers;
    variousIntegers.setSize(numSendBufs * 4);

    // We need various sets of integers to make MPI happy, each of which equal in
    // length to the number of processors:
    // - Send counts: how much data to send to each other process.
    // - Send cumsums: cumsums[i] indicates where to start reading from for the
    //   data to send. To processor `i`, we send the region starting at `sendBuf[cumsums[i]]`,
    //   extending for `counts[i]`-many elements.
    //
    // How much to send to each processor.
    Foam::UList<int> sendCounts{variousIntegers.begin(), numSendBufs};

    // Cumsum of sendCounts
    Foam::UList<int> sendCumsums{variousIntegers.begin() + numSendBufs, numSendBufs};

    // How much to receive from each processor.
    Foam::UList<int> recvCounts{variousIntegers.begin() + numSendBufs * 2, numSendBufs};

    // Cumsum of recvCounts.
    Foam::UList<int> recvCumsums{variousIntegers.begin() + numSendBufs * 3, numSendBufs};

    // Two loops to fill in the above in the obvious way...
    for (int i = 0; i < numSendBufs; i++) {
        label elementsToSend = sendBufs[i].size();
        int recordsToSend = (int) divRoundUp(elementsToSend, EltsPerRecord);

        sendCumsums[i] = numRecordsToSend;

        numRecordsToSend += recordsToSend;
        sendCounts[i] = recordsToSend;
    }

    for (int i = 0; i < numSendBufs; i++) {
        label elementsToRecv = recvSizes[i];
        int recordsToRecv = (int) divRoundUp(elementsToRecv, EltsPerRecord);

        recvBufs[i].setSize(elementsToRecv);

        recvCumsums[i] = numRecordsToRecv;
        numRecordsToRecv += recordsToRecv;
        recvCounts[i] = recordsToRecv;
    }

    // TODO: If we were smarter, we could avoid copying the data into these intermediate buffers.
    //       That'd require cooperation with the callsites, which might be optimistic.
    Foam::List<RecordTy> sendCoalesced;
    sendCoalesced.setSize(numRecordsToSend);
    Foam::List<RecordTy> recvCoalesced;
    recvCoalesced.setSize(numRecordsToRecv);

    // Combine the things to send into a single buffer.
    label offs = 0;
    for (int i = 0; i < numSendBufs; i++) {
        // `offs` counts in multiples of records, so the base pointer naturally ends up
        // record-aligned.
        // The copy-length is in terms of elements, naturally leaving the junk region at the end of
        // the last record untouched.
        const auto& buf = sendBufs[i];
        label numRecords = divRoundUp(buf.size(), EltsPerRecord);

        if
        (
            (buf.size() > 0)
         && (buf.cdata() != nullptr)
         && (sendCoalesced.data() != nullptr)
        )
        {
            memcpy(sendCoalesced.data() + offs, buf.cdata(), buf.size() * sizeof(T));
            offs += numRecords;
        }
    }

    // TODO: It seems that this routine is sometimes used even when what's happening is
    //       just master sending data to everyone. That could be more efficiently done
    //       with scatterv or such.
    // TODO: There's a suspiciously large number of calls to this function where the amount
    //       of data to both send and receive are zero! Depending how idiotically MPI
    //       handles that, there might be some savings to be had.
    MPI_Alltoallv(
        sendCoalesced.data(), sendCounts.data(), sendCumsums.data(), MpiRecordTy,
        recvCoalesced.data(), recvCounts.data(), recvCumsums.data(), MpiRecordTy,
        PstreamGlobals::MPICommunicators_[comm]
    );

    // Sprinkle the received values back out to the silly receive buffers.
    // We should refactor so the users just want a list of views into the coalesced buffers,
    // then we get to avoid some allocations/copies here.
    offs = 0;
    for (int i = 0; i < numSendBufs; i++) {
        label sz = recvSizes[i];
        label numRecords = divRoundUp(sz, EltsPerRecord);

        // As above, this memcpy just... ignores the junk zone.
        if
        (
            (sz > 0)
         && (recvBufs[i].data() != nullptr)
         && (recvCoalesced.data() != nullptr)
        )
        {
            memcpy(recvBufs[i].data(), recvCoalesced.data() + offs, sz * sizeof(T));
            offs += numRecords;
        }
    }
}

namespace Foam { // <- Only to make ADL work. It's basically an anonymous ns.

// TODO: exchangeBuf and exchangeContainer are essentially identical...
template<class T>
static void exchangeBuf
(
    const labelUList& sendSizes,
    const UList<const char*>& sendBufs,
    const labelUList& recvSizes,
    List<char*>& recvBufs,
    const label comm
)
{
    label startOfRequests = Pstream::nRequests();

    // Set up receives
    // ~~~~~~~~~~~~~~~
    forAll(recvSizes, proci)
    {
        if (proci != Pstream::myProcNo(comm) && recvSizes[proci] > 0)
        {
            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                recvBufs[proci],
                recvSizes[proci]*sizeof(T),
                0,
                comm
            );
        }
    }

    // Set up sends
    // ~~~~~~~~~~~~
    forAll(sendBufs, proci)
    {
        if (proci != Pstream::myProcNo(comm) && sendSizes[proci] > 0)
        {
            if
            (
               !UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    sendBufs[proci],
                    sendSizes[proci]*sizeof(T),
                    0,
                    comm
                )
            )
            {
                FatalErrorInFunction
                    << "Cannot send outgoing message. "
                    << "to:" << proci << " nBytes:"
                    << label(sendSizes[proci]*sizeof(T))
                    << Foam::abort(FatalError);
            }
        }
    }

    // Wait for all to finish
    // ~~~~~~~~~~~~~~~~~~~~~~
    Pstream::waitRequests(startOfRequests);
}

template<class Container, class T>
static void exchangeContainer
(
    const UList<Container>& sendBufs,
    const labelUList& recvSizes,
    List<Container>& recvBufs,
    const label comm
)
{
    label startOfRequests = Pstream::nRequests();

    // Set up receives
    // ~~~~~~~~~~~~~~~

    forAll(recvSizes, proci)
    {
        if (proci != Pstream::myProcNo(comm) && recvSizes[proci] > 0)
        {
            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                reinterpret_cast<char*>(recvBufs[proci].begin()),
                recvSizes[proci]*sizeof(T),
                0,
                comm
            );
        }
    }


    // Set up sends
    // ~~~~~~~~~~~~

    forAll(sendBufs, proci)
    {
        if (proci != Pstream::myProcNo(comm) && sendBufs[proci].size() > 0)
        {
            if
            (
               !UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    reinterpret_cast<const char*>(sendBufs[proci].begin()),
                    sendBufs[proci].size()*sizeof(T),
                    0,
                    comm
                )
            )
            {
                FatalErrorInFunction
                    << "Cannot send outgoing message. "
                    << "to:" << proci << " nBytes:"
                    << label(sendBufs[proci].size()*sizeof(T))
                    << Foam::abort(FatalError);
            }
        }
    }

    // Wait for all to finish
    // ~~~~~~~~~~~~~~~~~~~~~~
    Pstream::waitRequests(startOfRequests);
}

} // Namespace Foam

template<class Container, class T>
void Foam::Pstream::exchangePointToPoint
(
    const UList<Container>& sendBufs,
    const labelUList& recvSizes,
    List<Container>& recvBufs,
    label comm
)
{
    // Presize all receive buffers
    forAll(recvSizes, proci)
    {
        label nRecv = recvSizes[proci];
        if (proci != Pstream::myProcNo(comm) && nRecv > 0)
        {
            recvBufs[proci].setSize(nRecv);
        }
    }
    if (Pstream::maxCommsSize <= 0)
    {
        // Do the exchanging in one go
        exchangeContainer<Container, T>
        (
            sendBufs,
            recvSizes,
            recvBufs,
            comm
        );
    }
    else
    {
        // Determine the number of chunks to send. Note that we
        // only have to look at the sending data since we are
        // guaranteed that some processor's sending size is some other
        // processor's receive size. Also we can ignore any local comms.
        label maxNSend = 0;
        forAll(sendBufs, proci)
        {
            if (proci != Pstream::myProcNo(comm))
            {
                maxNSend = max(maxNSend, sendBufs[proci].size());
            }
        }
        const label maxNBytes = sizeof(T)*maxNSend;
        // We need to send maxNBytes bytes so the number of iterations:
        //  maxNBytes                           iterations
        //  ---------                           ----------
        //  0                                   0
        //  1..maxCommsSize                     1
        //  maxCommsSize+1..2*maxCommsSize      2
        //      etc.
        label nIter;
        if (maxNBytes == 0)
        {
            nIter = 0;
        }
        else
        {
            nIter = (maxNBytes-1)/Pstream::maxCommsSize+1;
        }
        reduce(nIter, maxOp<label>(), comm);
        List<const char*> charSendBufs(sendBufs.size());
        List<char*> charRecvBufs(sendBufs.size());
        labelList nRecv(sendBufs.size());
        labelList startRecv(sendBufs.size(), 0);
        labelList nSend(sendBufs.size());
        labelList startSend(sendBufs.size(), 0);
        for (label iter = 0; iter < nIter; iter++)
        {
            forAll(sendBufs, proci)
            {
                nSend[proci] = min
                (
                    Pstream::maxCommsSize,
                    sendBufs[proci].size()-startSend[proci]
                );
                charSendBufs[proci] =
                (
                    nSend[proci] > 0
                  ? reinterpret_cast<const char*>
                    (
                        &(sendBufs[proci][startSend[proci]])
                    )
                  : nullptr
                );
                nRecv[proci] = min
                (
                    Pstream::maxCommsSize,
                    recvBufs[proci].size()-startRecv[proci]
                );
                charRecvBufs[proci] =
                (
                    nRecv[proci] > 0
                  ? reinterpret_cast<char*>
                    (
                        &(recvBufs[proci][startRecv[proci]])
                    )
                  : nullptr
                );
            }
            exchangeBuf<T>
            (
                nSend,
                charSendBufs,
                nRecv,
                charRecvBufs,
                comm
            );
            forAll(nSend, proci)
            {
                startSend[proci] += nSend[proci];
                startRecv[proci] += nRecv[proci];
            }
        }
    }

    // Do myself
    recvBufs[Pstream::myProcNo(comm)] = sendBufs[Pstream::myProcNo(comm)];
}

template<class Container>
void Foam::Pstream::exchangeSizes
(
    const Container& sendBufs,
    labelList& recvSizes,
    label comm
)
{
    // sendBufs: what I want to send to everyone else.
    // recvSizes: populated with the amount everyone else wants to send to me.
    if (sendBufs.size() != UPstream::nProcs(comm))
    {
        FatalErrorInFunction
            << "Size of container " << sendBufs.size()
            << " does not equal the number of processors "
            << UPstream::nProcs(comm)
            << Foam::abort(FatalError);
    }

    labelList sendSizes(sendBufs.size());
    forAll(sendBufs, proci)
    {
        sendSizes[proci] = sendBufs[proci].size();
    }

    // When there's no parallelism, just copy the input to the output.
    if (UPstream::nProcs(comm) == 1 || !UPstream::parRun()) {
        recvSizes = std::move(sendSizes);
        return;
    }

    recvSizes.setSize(sendSizes.size());
    auto MpiRecordTy = UPstream::dataTypes[sizeof(label)];

    // Every process sends the amount it wants to send to every other process to that process.
    // `recvSizes` ends up containing the amount that process wants to send you.
    MPI_Alltoall(
        sendSizes.begin(), 1, MpiRecordTy,
        recvSizes.begin(), 1, MpiRecordTy,
        PstreamGlobals::MPICommunicators_[comm]
    );
}


template<class Container, class T>
void Foam::Pstream::exchange
(
    const UList<Container>& sendBufs,
    List<Container>& recvBufs,
    const label comm
)
{
    labelList recvSizes;
    exchangeSizes(sendBufs, recvSizes, comm);

    exchange<Container, T>(sendBufs, recvSizes, recvBufs, comm);
}


// ************************************************************************* //
