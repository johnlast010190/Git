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
    writeFieldsOftenFunctionObject

Description
    Lets the user specify an additional set of fields that is to be written
    more frequently than specified in the controlDict

    writeControl accepts the same values

Example usage: writes Temperature every third timestep after t=11

functions
(
    writeT {
        functionObjectLibs ("libsimpleFunctionObjects.so");
        type writeFieldsOften;
        fieldNames ( T );
        after 11;
        writeControl timeStep;
        writeIntervall 3;
    }
);

SourceFiles
    writeFieldsOftenFunctionObject.C

Contributors/Copyright:
    2008-2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2021 Engys Ltd

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef writeFieldsOftenFunctionObject_H
#define writeFieldsOftenFunctionObject_H

#include "writeFields/writeFieldsGeneral/writeFieldsGeneralFunctionObject.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class writeFieldsOftenFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class writeFieldsOftenFunctionObject
:
    public writeFieldsGeneralFunctionObject
{
protected:

    // Protected Data

        //- How writting is controlled
        Time::writeControls writeControl_;

        //- How often we write
        scalar writeInterval_;

        //- Store the last output
        label outputTimeIndex_;

        scalar baselineTime_;


    // Protected Member Functions

        //- Overwrite the default output time
        virtual bool outputTime(const bool forceWrite);


public:

    //- Runtime type information
    TypeName("writeFieldsOften");


    // Constructors

        //- Construct from components
        writeFieldsOftenFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        writeFieldsOftenFunctionObject
        (
            const writeFieldsOftenFunctionObject&
        ) = delete;


    // Member Functions

        bool start();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const writeFieldsOftenFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
