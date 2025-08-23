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
    writeAdditionalFieldsFunctionObject

Description
    Lets the user specify an additional set of fields that is to be written with the other
    fields the are written by AUTO_WRITE

Example usage: writes Enthalpy and viscosity after t=15

functions
(
    writeEnthalpy
    {
        functionObjectLibs ("libsimpleFunctionObjects.so");
        type writeAdditionalFields;
        fieldNames ( h mu );
        after 15;
    }
);

SourceFiles
    writeAdditionalFieldsFunctionObject.C

Contributors/Copyright:
    2008-2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef writeAdditionalFieldsFunctionObject_H
#define writeAdditionalFieldsFunctionObject_H

#include "writeFields/writeFieldsGeneral/writeFieldsGeneralFunctionObject.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class writeAdditionalFieldsFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class writeAdditionalFieldsFunctionObject
:
    public writeFieldsGeneralFunctionObject
{
public:

    //- Runtime type information
    TypeName("writeAdditionalFields");


    // Constructors

        //- Construct from components
        writeAdditionalFieldsFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        writeAdditionalFieldsFunctionObject
        (
            const writeAdditionalFieldsFunctionObject&
        ) = delete;


    // Member Functions

        bool start();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const writeAdditionalFieldsFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
