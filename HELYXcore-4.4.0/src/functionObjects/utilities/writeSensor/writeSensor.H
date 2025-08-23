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
    (c) 2019-2024 Engys Ltd.

Author
    2019. Nikolaos Magoulas (Engys Ltd.). All rights reserved.

Class
    Foam::functionObjects::writeSensor

Description
    Function object to write the fields of sensors.

Usage
    Example of function object specification:
    \verbatim
    writeSensor1
    {
        type        writeSensor;
        libs        ("libutilityFunctionObjects.so");

        sensor
        {
            type    magGradU;
            writeSensorField    true; // Don't forget this setting
            ...
        }
    }
    \endverbatim

SourceFiles
    writeSensor.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writeSensor_H
#define functionObjects_writeSensor_H

#include "db/functionObjects/functionObject/functionObject.H"
#include "cfdTools/general/sensor/sensor/sensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class writeSensor Declaration
\*---------------------------------------------------------------------------*/

class writeSensor
:
    public functionObject
{
private:

    // Private data

        //- Reference to Db
        const objectRegistry& obr_;

        //- Sensor function of whose the sensorField is written
        autoPtr<sensor<vector>> sensor_;


public:

    //- Runtime type information
    TypeName("writeSensor");


    // Constructors

        //- Construct from Time and dictionary
        writeSensor
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        writeSensor(const writeSensor&) = delete;


    //- Destructor
    virtual ~writeSensor();


    // Member Functions

        //- Read the writeSensor data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Write the registered objects
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const writeSensor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
