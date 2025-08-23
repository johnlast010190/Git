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
    (c) 1991-2008 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::shiftFieldConstantVectorPluginFunction

Description

SourceFiles
    shiftFieldConstantVectorPluginFunction.C

Contributors/Copyright:
    2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef shiftFieldConstantVectorPluginFunction_H
#define shiftFieldConstantVectorPluginFunction_H

#include "shiftFieldGeneralPluginFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class shiftFieldConstantVectorPluginFunction Declaration
\*---------------------------------------------------------------------------*/

template<class Type,meshToMeshOrder Order>
class shiftFieldConstantVectorPluginFunction
:
    public shiftFieldGeneralPluginFunction<Type, Order>
{
    // Private Data

        vector shiftValue_;


protected:

    // Protected Member Functions

        tmp<pointField> displacement();

        void setArgument(label index, const vector& content);


public:

    //- Runtime type information
    TypeName("shiftFieldConstantVectorPluginFunction");


    // Constructors

        //- Construct from components
        shiftFieldConstantVectorPluginFunction
        (
            const FieldValueExpressionDriver& parentDriver,
            const word& name
        );

        //- Disallow default bitwise copy construction
        shiftFieldConstantVectorPluginFunction
        (
            const shiftFieldConstantVectorPluginFunction&
        ) = delete;


    //- Destructor
    virtual ~shiftFieldConstantVectorPluginFunction()
    {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const shiftFieldConstantVectorPluginFunction&) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
