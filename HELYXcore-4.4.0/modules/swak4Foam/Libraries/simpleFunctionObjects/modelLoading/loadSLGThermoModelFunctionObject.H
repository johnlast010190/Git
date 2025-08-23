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
    loadSLGThermoModelFunctionObject

Description
    Loads a SLGThermo model into memory

SourceFiles
    loadSLGThermoModelFunctionObject.C

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef loadSLGThermoModelFunctionObject_H
#define loadSLGThermoModelFunctionObject_H

#include "modelLoadingFunctionObject.H"
#include "SLGThermo/SLGThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class loadSLGThermoModelFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class loadSLGThermoModelFunctionObject
:
    public modelLoadingFunctionObject<SLGThermo>
{
protected:

    // Protected Member Functions

        virtual autoPtr<SLGThermo> initModel();


public:

    //- Runtime type information
    TypeName("loadSLGThermoModelFunctionObject");


    // Constructors

        //- Construct from components
        loadSLGThermoModelFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        loadSLGThermoModelFunctionObject
        (
            const loadSLGThermoModelFunctionObject&
        ) = delete;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const loadSLGThermoModelFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
