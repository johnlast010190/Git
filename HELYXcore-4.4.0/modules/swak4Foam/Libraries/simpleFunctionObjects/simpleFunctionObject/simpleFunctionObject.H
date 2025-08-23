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
    simpleFunctionObject

Description
    Basis for the other function objects here. Evaluates the most common options

SourceFiles
    simpleFunctionObject.C

Contributors/Copyright:
    2008-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef simpleFunctionObject_H
#define simpleFunctionObject_H

#include "db/functionObjects/functionObject/functionObject.H"
#include "db/dictionary/dictionary.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"
#include "primitives/bools/lists/boolList.H"
#include "fields/Fields/Field/Field.H"
#include "primitives/enums/NamedEnum.H"
#include "include/swak.H"
#include "helpers/switchableNotImplemented.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class simpleFunctionObject Declaration
\*---------------------------------------------------------------------------*/

    class objectRegistry;

class simpleFunctionObject
:
    public functionObject
{
public:

    enum outputControlModeType
    {
        ocmTimestep,
        ocmDeltaT,
        ocmOutputTime,
        ocmStartup,
        ocmOutputTimeAndStartup,
        ocmTimestepAndStartup,
        ocmDeltaTAndStartup
    };


private:

    // Private data

        //- Do output on the terminal
        bool verbose_;

        //- Write debug-information to the terminal
        bool writeDebug_;

        //- Time after which calculation should start
        scalar after_;

        //- How many time-steps since last output
        label timeSteps_;

        //- Possible outputControl-modes
        static const NamedEnum<outputControlModeType,7> outputControlModeTypeNames_;
        outputControlModeType outputControlMode_;

        //- How many time-steps between outputs
        label outputInterval_;

        //- DeltaT for output
        scalar outputDeltaT_;

        const Time& time_;

        scalar lastWrite_;

        bool started_;

#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
        //- Make sure that execute is not called twice at one time-step
        label lastTimeStepExecute_;

        //- Allow executing more than once - for instance in fvOptions
        bool executeMoreThanOnce_;


    // Private Member Functions

        bool ensureExecuteOnce();
#endif

protected:

    // Protected Data

        dictionary dict_;

        word regionName_;

        string regionString_;

        const objectRegistry& obr_;


    // Protected Member Functions

        outputControlModeType outputControlMode() const
        {
            return outputControlMode_;
        }

        const char* regionString() const
        {
            return regionString_.c_str();
        }

        //- Do output on the terminal
        bool verbose() const
        {
            return verbose_;
        }

        //- Write debug-information to tthe terminal
        bool writeDebug() const
        {
            return writeDebug_;
        }

        //- Do output this time?
        virtual bool outputTime(const bool forceWrite);


public:

    //- Runtime type information
    TypeName("simpleFunctionObject");


    // Constructors

        //- Construct from components
        simpleFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        simpleFunctionObject(const simpleFunctionObject&) = delete;


    // Member Functions

        //- Start is called at the start of the time-loop
        virtual bool start();

        //- Execute is called at each ++ or += of the time-loop
        virtual bool execute(const bool forceWrite);

#ifdef FOAM_FUNCTIONOBJECT_EXECUTE_HAS_NO_FORCE
        bool execute()
        {
            return execute(false);
        }
#endif

#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
        //- Interface for new functionObject
        virtual bool execute();
        virtual bool write();
#endif

        //- Read and set the function object if it's data has changed
        virtual bool read(const dictionary& dict);

        //- Sample and write
        virtual void writeSimple() = 0;

#ifdef FOAM_FUNCTIONOBJECT_NEEDS_UPDATEMESH
        //- Make functionObject happy but fail
        //  because this is not an implementation
        virtual void topoChange(const polyTopoChangeMap&)
        {
            switchableNotImplemented
            (
                "simpleFunctionObject::topoChange",
                dict()
            );
        }

        //- Make functionObject happy but fail
        //  because this is not an implementation
        virtual void movePoints(const polyMesh&)
        {
            switchableNotImplemented
            (
                "simpleFunctionObject::movePoints",
                dict()
            );
        }
#endif

        //- Optional stuff to do after write
        virtual void flush();

        //- the actual time
        const Time& time() const
        {
            return time_;
        }

        //- The actual time
        scalar after() const
        {
            return after_;
        }

        const objectRegistry& obr() const
        {
            return obr_;
        }

        const dictionary& dict() const
        {
            return dict_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const simpleFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
