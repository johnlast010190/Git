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
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "rtppFunctionObject.H"

#include "signals/sigFpe.H"

#include "PostDictObjectProviderDatabase/postDictObjectProviderDatabase.H"
#include "PostDictObjectProviderDatabase/postDict/postDictKeys.H"
#include "Utils/memoryMonitor.H"

#if BASIC_PROFILE == 1
#include "Utils/basicProfiler.H"
static int anExecutions = 0;
#endif

#undef Log
#include "vtkLogger.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::rtppFunctionObject::rtppFunctionObject
    (
        const word &name,
        const Time &runTime,
        const dictionary& dict
    )
    :
    stateFunctionObject(name, runTime)
{
#if BASIC_PROFILE == 1
    runTimeVis::basicProfiler::initProfiler();
    runTimeVis::basicProfiler::startGroup("creating function object");
#endif

    rtppFunctionObject::read(dict);

    Foam::functionObjects::runTimeVis::PostDictObjectProviderDatabase::initializeStaticCommon(debug);

    if (debug) vtkLogF(INFO, "Initializing %s. Memory consumption: current: %zu mb; peak: %zu mb", name.c_str(), MemoryMonitor::getCurrentRSSMByte(), MemoryMonitor::getPeakRSSMByte());
}

bool Foam::functionObjects::rtppFunctionObject::read(const dictionary& dict)
{
    bool returnValue = stateFunctionObject::read(dict);

    if (readDebugFlag(dict))
    {
        debug = 1;
        DebugInformation << "    debugging on" << endl;
    }

    return returnValue;
}

bool Foam::functionObjects::rtppFunctionObject::readDebugFlag(const Foam::dictionary &dictionary)
{
    return dictionary.lookupOrDefault(runTimeVis::controlDictKeys::DEBUG_KEY, false);
}

void Foam::functionObjects::rtppFunctionObject::finishInitialization() const
{
    if (debug) vtkLogF(INFO, "Finished initializing. Memory consumption: current: %zu mb; peak: %zu mb", MemoryMonitor::getCurrentRSSMByte(), MemoryMonitor::getPeakRSSMByte());
}

Foam::functionObjects::rtppFunctionObject::~rtppFunctionObject()
{
    Foam::functionObjects::runTimeVis::PostDictObjectProviderDatabase::finalizeStaticCommon();
};

// ************************************************************************* //
