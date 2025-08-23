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
    (c) 2024 Engys Ltd.

Class
    writeOldTimesOnSignalFunctionObject

Description
    This function object tries to record old fields and writes them in
    case of a Signal

SourceFiles
    writeOldTimesOnSignalFunctionObject.C

Contributors/Copyright:
    2014-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef writeOldTimesOnSignalFunctionObject_H
#define writeOldTimesOnSignalFunctionObject_H

#include "db/functionObjects/functionObject/functionObject.H"
#include "db/dictionary/dictionary.H"
#include "DataStructures/TimeClone/TimeCloneList.H"
#include "containers/Lists/DynamicList/DynamicList.H"
#include <signal.h>
#include "include/swak.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class writeOldTimesOnSignalFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class objectRegistry;

class writeOldTimesOnSignalFunctionObject
:
    public functionObject
{
    // Private Data

        TimeCloneList times_;
        const Time& theTime_;

    class SignalHandlerInfo
    {
    protected:

        // Protected Data

            struct sigaction oldAction_;
            word name_;
            int sig_;
            bool set_;


    public:

        // Constructors

            //- Construct from components
            SignalHandlerInfo(word name, int sig);

            //- Construct from components
            SignalHandlerInfo()
            :
                name_("nix"),
                sig_(-1),
                set_(false)
            {}

        //- Destructor
        ~SignalHandlerInfo()
        {};


        // Member Functions

            void resetHandler();

            int sig() const
            {
                return sig_;
            }

            bool set() const
            {
                return set_;
            }
    };

    // Privata Member Functions

        static void sigHandler(int);


    // Private Data

        static writeOldTimesOnSignalFunctionObject *singleton_;

        DynamicList<SignalHandlerInfo> handlers_;

        bool writeCurrent_;

        scalar sleepSecondsBeforeReraising_;

        bool sigFPE_;
        bool sigSEGV_;
        bool sigINT_;
        bool sigTERM_;
        bool sigQUIT_;
        bool sigUSR1_;
        bool sigUSR2_;

        bool alreadyDumped_;

        bool itWasMeWhoReraised_;


public:

    //- Runtime type information
    TypeName("writeOldTimesOnSignal");


    // Constructors

        //- Construct from components
        writeOldTimesOnSignalFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        writeOldTimesOnSignalFunctionObject
        (
            const writeOldTimesOnSignalFunctionObject&
        ) = delete;


    // Member Functions

        //- Start is called at the start of the time-loop
        virtual bool start();

        //- Execute is called at each ++ or += of the time-loop
#ifdef FOAM_FUNCTIONOBJECT_EXECUTE_HAS_NO_FORCE
        virtual bool execute();
#else
        virtual bool execute(const bool forceWrite);

#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
        virtual bool execute()
        {
            return execute(false);
        }

        virtual bool write()
        {
            return execute(true);
        }
#endif
#endif


        //- Read and set the function object if it's data has changed
        virtual bool read(const dictionary& dict);

#ifdef FOAM_FUNCTIONOBJECT_NEEDS_UPDATEMESH

        virtual void topoChange(const polyTopoChangeMap&)
        {}

        virtual void movePoints(const polyMesh&)
        {}
#endif


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const writeOldTimesOnSignalFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
