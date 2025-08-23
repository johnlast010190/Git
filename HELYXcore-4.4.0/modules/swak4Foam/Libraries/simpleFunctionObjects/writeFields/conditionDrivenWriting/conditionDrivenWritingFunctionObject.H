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
    conditionDrivenWritingFunctionObject

Description
    General class for functionObjects that write timesteps if certain conditions
    are met

SourceFiles
    conditionDrivenWritingFunctionObject.C

Contributors/Copyright:
    2008-2011, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id: conditionDrivenWritingFunctionObject.H,v 78b0d113b99b 2013-02-25 16:12:41Z bgschaid $
\*---------------------------------------------------------------------------*/

#ifndef conditionDrivenWritingFunctionObject_H
#define conditionDrivenWritingFunctionObject_H

#include "simpleFunctionObject/simpleFunctionObject.H"
#include "primitives/enums/NamedEnum.H"
#include "DataStructures/TimeClone/TimeCloneList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class conditionDrivenWritingFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class conditionDrivenWritingFunctionObject
:
    public simpleFunctionObject
{
public:

    enum writeControlModeType
    {
        scmWriteAlways,
        scmWriteNTimesteps,
        scmWriteIntervall,
        scmWriteUntilSwitch
    };

    enum cooldownModeType
    {
        cdmNoCooldown,
        cdmNTimesteps,
        cdmIntervall,
        cdmRetrigger
    };


private:

    // Private Data

        enum theStateType
        {
            stateWaiting,
            stateWriting,
            stateStartCooldown,
            stateCooldown
        };

        static const NamedEnum<writeControlModeType,4> writeControlModeTypeNames_;
        writeControlModeType writeControlMode_;

        static const NamedEnum<cooldownModeType,4> cooldownModeTypeNames_;
        cooldownModeType cooldownMode_;

        bool storeAndWritePreviousState_;

        theStateType theState_;

        label writeTimesteps_;
        scalar writeIntervall_;
        label cooldownTimesteps_;
        scalar cooldownIntervall_;

        label timestepForStateChange_;
        scalar timeForStateChange_;

        autoPtr<TimeCloneList> lastTimes_;


    // Private Member Functions

        void storePreviousState();

        void writePreviousState();

        void writeNow();

        bool alreadyWritten(word timename);

        bool checkWrite();

        bool checkCooldown();


protected:

    // Protected Member Functions

        //- should writing start now?
        virtual bool checkStartWriting() = 0;

        //- switch off writing? Does not have to be implemented
        virtual bool checkStopWriting();

        //- switch off cooldown? Does not have to be implemented
        virtual bool checkStopCooldown();

        cooldownModeType cooldownMode()
        {
            return cooldownMode_;
        }

        writeControlModeType writeControlMode()
        {
            return writeControlMode_;
        }


public:

    //- Runtime type information
    TypeName("conditionDrivenWriting");


    // Constructors

        //- Construct from components
        conditionDrivenWritingFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        conditionDrivenWritingFunctionObject
        (
            const conditionDrivenWritingFunctionObject&
        ) = delete;


    // Member Functions

        bool start();

        void writeSimple();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const conditionDrivenWritingFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
