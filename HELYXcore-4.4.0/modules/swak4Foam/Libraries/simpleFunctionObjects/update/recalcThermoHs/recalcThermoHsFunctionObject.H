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
    recalcThermoHsFunctionObject

Description
    Recalculate sensible enthalpy according to the current temperature

SourceFiles
    recalcThermoHsFunctionObject.C

Contributors/Copyright:
    2008-2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef recalcThermoHsFunctionObject_H
#define recalcThermoHsFunctionObject_H

#include "update/updateSimple/updateSimpleFunctionObject.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class recalcThermoHsFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class recalcThermoHsFunctionObject
:
    public updateSimpleFunctionObject
{
protected:

    // Protected Member Functions

        //- Do the actual calculation
        void recalc();


public:

    //- Runtime type information
    TypeName("recalcThermoHs");

    // Constructors

        //- Construct from components
        recalcThermoHsFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        recalcThermoHsFunctionObject
        (
            const recalcThermoHsFunctionObject&
        ) = delete;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const recalcThermoHsFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
