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
    readGravitationFunctionObject

Description
    Read the gravitational acceleration. This is done to satisfy boundary
    conditions that rely on it for solvers/utilities that don't provide it

SourceFiles
    readGravitationFunctionObject.C

Contributors/Copyright:
    2008-2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef readGravitationFunctionObject_H
#define readGravitationFunctionObject_H

#include "simpleFunctionObject/simpleFunctionObject.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class readGravitationFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class readGravitationFunctionObject
:
    public simpleFunctionObject
{
    // Private Data

        //- The actual field
        uniformDimensionedVectorField g_;


public:

    //- Runtime type information
    TypeName("readGravitation");


    // Constructors

        //- Construct from components
        readGravitationFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        readGravitationFunctionObject
        (
            const readGravitationFunctionObject&
        ) = delete;


    // Member Functions

        void writeSimple();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const readGravitationFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
