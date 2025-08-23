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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::outputFilterOutputControl

Description
    An output control for function objects.
    The default is time-step execution at every interval.

SourceFiles
    outputFilterOutputControl.C

\*---------------------------------------------------------------------------*/

#ifndef outputFilterOutputControl_H
#define outputFilterOutputControl_H

#include "db/dictionary/dictionary.H"
#include "db/Time/Time.H"
#include "primitives/enums/NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class outputFilterOutputControl Declaration
\*---------------------------------------------------------------------------*/

class outputFilterOutputControl
{
public:

    //- The output control options
    enum outputControls
    {
        ocTimeStep,       /*!< execution is coupled to the time-step */
        ocOutputTime,     /*!< execution is coupled to the output-time */
        ocAdjustableTime, /*!< Adjust time step for dumping */
        ocRunTime,        /*!< run time for dumping */
        ocClockTime,      /*!< clock time for dumping */
        ocCpuTime,        /*!< cpu time for dumping */
        ocNone            /*!< no output */
    };


private:

    // Private data

        //- Time object
        const Time& time_;

        //- Prefix
        const word prefix_;

        //- String representation of outputControls enums
        static const NamedEnum<outputControls, 7> outputControlNames_;

        //- Type of output
        outputControls outputControl_;

        //- The execution interval (in time steps) when using \c timeStep mode,
        //  a value <= 1 means execute at every time step
        label outputInterval_;

        //- Dumping counter for ocOutputTime or index dump for ocAdjustableTime
        label outputTimeLastDump_;

        //- Dump each deltaT (adjust Ttime)
        scalar writeInterval_;


public:

    // Constructors

        //- Construct from Time object and dictionary
        outputFilterOutputControl
        (
            const Time&,
            const dictionary&,
            const word& prefix
        );

        //- Disallow default bitwise copy construction
        outputFilterOutputControl(const outputFilterOutputControl&) = delete;


    //- Destructor
    ~outputFilterOutputControl();


    // Member Functions

        //- Read from dictionary
        void read(const dictionary&);

        //- Return const access to the Time object
        const Time& time() const
        {
            return time_;
        }

        //- Flag to indicate whether to output
        bool output();

        //- Return outputControl
        outputControls outputControl() const
        {
            return outputControl_;
        }

        //- Return writeInterval
        scalar writeInterval() const
        {
            return writeInterval_;
        }

        //- Return outputTimeLastDump
        label outputTimeLastDump() const
        {
            return outputTimeLastDump_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const outputFilterOutputControl&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
