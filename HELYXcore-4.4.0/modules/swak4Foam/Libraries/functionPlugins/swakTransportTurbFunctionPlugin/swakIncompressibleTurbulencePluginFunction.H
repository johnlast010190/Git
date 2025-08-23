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
    (c) 1991-2008 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::swakIncompressibleTurbulencePluginFunction

Description

SourceFiles
    swakIncompressibleTurbulencePluginFunction.C

Contributors/Copyright:
    2012-2013, 2015-2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef swakIncompressibleTurbulencePluginFunction_H
#define swakIncompressibleTurbulencePluginFunction_H

#include "swakTransportModelsPluginFunction.H"

#include "include/swakTurbulence.H"

#ifdef FOAM_NEW_TURBULENCE_STRUCTURE
#include "singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#else
#include "turbulenceModel.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class swakIncompressibleTurbulencePluginFunction Declaration
\*---------------------------------------------------------------------------*/
class swakIncompressibleTurbulencePluginFunction
:
    public swakTransportModelsPluginFunction
{
protected:

    // Protected Member Functions

        static const incompressible::turbulenceModel& turbInternal
        (
            const fvMesh& reg
        );

        const incompressible::turbulenceModel& turb();


public:

    //- Runtime type information
    TypeName("swakIncompressibleTurbulencePluginFunction");

    // Constructors

        //- Construct from components
        swakIncompressibleTurbulencePluginFunction(
            const FieldValueExpressionDriver &parentDriver,
            const word &name,
            const word &returnValueType
        );

        //- Disallow default bitwise copy construction
        swakIncompressibleTurbulencePluginFunction
        (
            const swakIncompressibleTurbulencePluginFunction&
        ) = delete;


    //- Destructor
    virtual ~swakIncompressibleTurbulencePluginFunction()
    {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const swakIncompressibleTurbulencePluginFunction&
        ) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
