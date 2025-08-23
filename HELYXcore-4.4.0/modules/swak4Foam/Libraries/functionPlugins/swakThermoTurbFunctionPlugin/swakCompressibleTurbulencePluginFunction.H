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
    Foam::swakCompressibleTurbulencePluginFunction

Description

SourceFiles
    swakCompressibleTurbulencePluginFunction.C

Contributors/Copyright:
    2012-2013, 2015-2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef swakCompressibleTurbulencePluginFunction_H
#define swakCompressibleTurbulencePluginFunction_H

#include "swakThermophysicalPluginFunction.H"

#include "include/swakTurbulence.H"

#ifdef FOAM_NEW_TURBULENCE_STRUCTURE
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#else
#include "turbulenceModel.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class swakCompressibleTurbulencePluginFunction Declaration
\*---------------------------------------------------------------------------*/
class swakCompressibleTurbulencePluginFunction
:
    public swakThermophysicalPluginFunction<swakFluidThermoType>
{
protected:

    // Protected Member Functions

        static const compressible::turbulenceModel& turbInternal
        (
            const fvMesh& reg
        );

        const compressible::turbulenceModel& turb();


public:

    //- Runtime type information
    TypeName("swakCompressibleTurbulencePluginFunction");


    // Constructors

        //- Construct from components
        swakCompressibleTurbulencePluginFunction
        (
            const FieldValueExpressionDriver& parentDriver,
            const word& name,
            const word& returnValueType
        );

        //- Disallow default bitwise copy construction
        swakCompressibleTurbulencePluginFunction
        (
            const swakCompressibleTurbulencePluginFunction&
        ) = delete;


    //- Destructor
    virtual ~swakCompressibleTurbulencePluginFunction() {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const swakCompressibleTurbulencePluginFunction&
        ) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
