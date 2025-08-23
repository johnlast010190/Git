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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "ParallelUtils.H"

#include "vtkMultiProcessController.h"
#include "vtkCellData.h"
#include "vtkImplicitFunction.h"
#include "vtkDoubleArray.h"

#include "db/IOstreams/Pstreams/Pstream.H"

namespace Foam::functionObjects::runTimeVis::ParallelUtils
{
void tradeValueWithProc(const scalar *input, scalar *output, label length, label procId)
{
    if (isRunningInParallel())
    {
        vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();
        int thisProc = controller->GetLocalProcessId();

        if (procId > thisProc)
        {
            controller->Receive(output, length, procId, 43);
            controller->Send(input, length, procId, 44);
        }
        else
        {
            controller->Send(input, length, procId, 43);
            controller->Receive(output, length, procId, 44);
        }
    }
    else
    {
        std::memcpy(output, input, length);
    }
}

bool isRunningInParallel()
{
    return Pstream::parRun();
}

bool isMaster()
{
    return Pstream::master();
}

int localProcessId()
{
    return Pstream::myProcNo();
}

label getMaxLabelFromAllProcs(label value, vtkMultiProcessController *controller)
{
    if (controller)
    {
        int np = controller->GetNumberOfProcesses();
        std::vector<label> values(np);
        controller->AllGather(&value, values.data(), 1);
        for (int i = 0; i < np; i++)
        {
            if (value < values[i])
            {
                value = values[i];
            }
        }
    }
    return value;
}

void debugInfoParallel(const std::string &s)
{
    int procId = 0;
    if (isRunningInParallel())
    {
        vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();
        procId = controller->GetLocalProcessId();
    }
    printf("Debug [%i]: %s \n", procId, s.c_str());
}

void debugInfoAllProcs(const string &text, scalar value)
{
    if (isRunningInParallel())
    {
        vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();
        int nProcs = controller->GetNumberOfProcesses();

        std::vector<scalar> allScalars(nProcs);

        controller->Gather(&value, allScalars.data(), 1, 0);

        if (controller->GetLocalProcessId() == 0)
        {
            Foam::Info<< text;
            for (int i = 0; i < nProcs; i++)
            {
                Foam::Info << ";Proc " << i << ": " << allScalars[i];
            }
            Foam::Info<< Foam::endl;
        }
    }
    else
    {
        Foam::Info<< text << value << endl;
    }
}

void debugInfoAllProcs(const string &text, const string &value)
{
    if (isRunningInParallel())
    {
        vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();
        int nProcs = controller->GetNumberOfProcesses();
        int localProcessId = controller->GetLocalProcessId();

        if (localProcessId == 0)
        {
            Foam::Info<< text;
            Foam::Info<< ";Proc 0: " + value;
            int stringLength = 0;
            for (int i = 1; i < nProcs; i++)
            {
                controller->Receive(&stringLength, 1, i, 6743);
                std::vector<char> receivedString(stringLength + 1);
                memset(receivedString.data(), 0, stringLength + 1);
                controller->Receive(receivedString.data(), stringLength, i, 2326);
                Foam::Info<< ";Proc " << i << ": " + string(receivedString.data());
            }
            Foam::Info<< Foam::endl;
        }
        else
        {
            auto stringLength = static_cast<vtkIdType>(value.length());
            controller->Send(&stringLength, 1, 0, 6743);
            controller->Send(value.c_str(), stringLength, 0, 2326);
        }
    }
    else
    {
        Foam::Info<< text << value << endl;
    }
}

List<List<scalar>> allGatherScalarList(const List<scalar> &values)
{
    debugInfoParallel("gathering");
    List<List<scalar>> listList = gatherScalarList(values);
    debugInfoParallel("scattering");
    listList = scatterScalarListList(listList);
    debugInfoParallel("done");
    return listList;
}

List<List<scalar>> gatherScalarList(const List<scalar> &values)
{
    if (isRunningInParallel())
    {
        vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();
        int nProcs = controller->GetNumberOfProcesses();
        int thisProc = controller->GetLocalProcessId();

        int length = values.size();
        if (thisProc == 0)
        {
            List<List<scalar>> procPointDistances(nProcs);
            procPointDistances[0] = values;
            for (int i = 1; i < nProcs; i++)
            {
                controller->Receive(&length, 1, i, 1646);
                List<scalar> received(length);
                controller->Receive(received.data(), length, i, 4879);
                procPointDistances[i] = received;
            }
            return procPointDistances;
        }
        else
        {
            controller->Send(&length, 1, 0, 1646);
            controller->Send(values.cdata(), length, 0, 4879);
            return {};
        }
    }
    else
    {
        List<List<scalar>> procPointDistances(1);
        procPointDistances[0] = values;
        return procPointDistances;
    }
}

List<List<scalar>> scatterScalarListList(const List<List<scalar>> &values)
{
    if (isRunningInParallel())
    {
        vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();
        int nProcs = controller->GetNumberOfProcesses();
        int thisProc = controller->GetLocalProcessId();

        if (thisProc == 0)
        {
            int nLists = values.size();
            for (int remoteProcessId = 1; remoteProcessId < nProcs; remoteProcessId++)
            {
                controller->Send(&nLists, 1, remoteProcessId, 1647);
                for (int listId = 0; listId < nLists; listId++)
                {
                    int length = values[listId].size();
                    controller->Send(&length, 1, remoteProcessId, 1648);
                    controller->Send(values[listId].cdata(), length, remoteProcessId, 4880);
                }
            }
            return values;
        }
        else
        {
            int nLists;
            controller->Receive(&nLists, 1, 0, 1647);
            List<List<scalar>> receivedLists(nLists);
            for (int listId = 0; listId < nLists; listId++)
            {
                int length;
                controller->Receive(&length, 1, 0, 1648);
                List<scalar> receivedList(length);
                controller->Receive(receivedList.data(), length, 0, 4880);
                receivedLists[listId] = receivedList;
            }
            return receivedLists;
        }
    }
    else
    {
        return values;
    }
}

}
