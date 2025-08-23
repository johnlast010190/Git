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
    loadCompressibleTurbulenceModelFunctionObject

Description
    Loads a turbulence model into memory

SourceFiles
    loadCompressibleTurbulenceModelFunctionObject.C

Contributors/Copyright:
    2012-2013, 2015-2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef loadCompressibleTurbulenceModelFunctionObject_H
#define loadCompressibleTurbulenceModelFunctionObject_H

#include "modelLoadingFunctionObject.H"

#include "include/swakTurbulence.H"

#ifdef FOAM_NEW_TURBULENCE_STRUCTURE
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#else
#include "compressible/turbulenceModel/turbulenceModel.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class loadCompressibleTurbulenceModelFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class loadCompressibleTurbulenceModelFunctionObject
:
    public modelLoadingFunctionObject<compressible::turbulenceModel>
{
protected:

    // Protected Member Functions

        virtual autoPtr<compressible::turbulenceModel> initModel();


public:

    //- Runtime type information
    TypeName("loadCompressibleTurbulenceModelFunctionObject");


    // Constructors

        //- Construct from components
        loadCompressibleTurbulenceModelFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        loadCompressibleTurbulenceModelFunctionObject
        (
            const loadCompressibleTurbulenceModelFunctionObject&
        ) = delete;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const loadCompressibleTurbulenceModelFunctionObject&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
