/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) ICE Stroemungsfoschungs GmbH

Contributors/Copyright:
    2014-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "provokeSignalFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "include/swakTime.H"

#include <signal.h>

#include "containers/HashTables/HashSet/HashSet.H"

#include "db/objectRegistry/objectRegistry.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(provokeSignalFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        provokeSignalFunctionObject,
        dictionary
    );

    template<>
    const char* NamedEnum<Foam::provokeSignalFunctionObject::possibleSignals,7>::names[]=
    {
        "FPE",
        "SEGV",
        "INT",
        "TERM",
        "QUIT",
        "USR1",
        "USR2"
    };
    const NamedEnum<provokeSignalFunctionObject::possibleSignals,7> provokeSignalFunctionObject::possibleSignalsNames_;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

provokeSignalFunctionObject::provokeSignalFunctionObject
(
    const word &name,
    const Time& t,
    const dictionary& dict
)
:
    simpleFunctionObject(name,t,dict),
    signalToRaise_(
        possibleSignalsNames_[
            word(dict.lookup("signalToRaise"))
        ]
    ),
    timeToRaise_(
        dict.lookup<scalar>("timeToRaise")
    ),
    raiseOnThisProc_(false)
{
    if (Pstream::parRun()) {
        HashSet<label> allProcs(
            labelList(
                dict.lookup("processorsToRaiseSignal")
            )
        );
        raiseOnThisProc_=allProcs.found(
            Pstream::myProcNo()
        );
    } else {
        if (!dict.found("processorsToRaiseSignal")) {
            WarningIn("provokeSignalFunctionObject::provokeSignalFunctionObject")
                << "No entry 'processorsToRaiseSignal' in " << dict.name()
                << nl << "Not needed now but needed in parallel runs" << endl;
        }
        raiseOnThisProc_=true;
    }

    if (raiseOnThisProc_) {
        Pout<< endl;
        Pout<< "Will raise signal " << possibleSignalsNames_[signalToRaise_]
            << " at time " << timeToRaise_
            << " and there is nothing you can do about it. "
            << "In fact the only possible use of this is testing" << endl;
        Pout<< endl;
    }
}

void provokeSignalFunctionObject::writeSimple()
{
    if (
        raiseOnThisProc_
        &&
        obr_.time().value()>=timeToRaise_
    ) {
        Pout<< endl
            << "The time has come. Raising "
            << possibleSignalsNames_[signalToRaise_] << endl
            << "It was nice knowing you"
            << endl;


       switch(signalToRaise_) {
            case sigFPE:
                raise(SIGFPE);
                break;
            case sigSEGV:
                raise(SIGSEGV);
                break;
            case sigINT:
                raise(SIGINT);
                break;
            case sigTERM:
                raise(SIGTERM);
                break;
            case sigQUIT:
                raise(SIGQUIT);
                break;
            case sigUSR1:
                raise(SIGUSR1);
                break;
            case sigUSR2:
                raise(SIGUSR2);
                break;
            default:
                FatalErrorIn("provokeSignalFunctionObject::writeSimple()")
                    << "Unimplemented signal "
                        << possibleSignalsNames_[signalToRaise_]
                        << endl;
        }
        Pout<< "We signaled (this should not be seen)" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

} // namespace Foam

// ************************************************************************* //
