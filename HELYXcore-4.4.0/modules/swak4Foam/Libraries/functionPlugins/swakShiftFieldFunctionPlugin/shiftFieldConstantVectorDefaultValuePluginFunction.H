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
    Foam::shiftFieldConstantVectorDefaultValuePluginFunction

Description

SourceFiles
    shiftFieldConstantVectorDefaultValuePluginFunction.C

Contributors/Copyright:
    2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef shiftFieldConstantVectorDefaultValuePluginFunction_H
#define shiftFieldConstantVectorDefaultValuePluginFunction_H

#include "shiftFieldGeneralPluginFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class shiftFieldConstantVectorDefaultValuePluginFunction Declaration
\*---------------------------------------------------------------------------*/

template<class Type,meshToMeshOrder Order>
class shiftFieldConstantVectorDefaultValuePluginFunction
:
    public shiftFieldGeneralPluginFunction<Type, Order>
{
    // Private Data

        vector shiftValue_;


protected:

    // Protected Member Data

        tmp<pointField> displacement();

        void setArgument(label index, const Type& content);

        void setArgument(label index, const vector& content);


public:

    //- Runtime type information
    TypeName("shiftFieldConstantVectorDefaultValuePluginFunction");

    // Constructors

        //- Construct from components
        shiftFieldConstantVectorDefaultValuePluginFunction
        (
            const FieldValueExpressionDriver& parentDriver,
            const word& name
        );

        //- Disallow default bitwise copy construction
        shiftFieldConstantVectorDefaultValuePluginFunction
        (
            const shiftFieldConstantVectorDefaultValuePluginFunction&
        ) = delete;


    //- Destructor
    virtual ~shiftFieldConstantVectorDefaultValuePluginFunction()
    {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const shiftFieldConstantVectorDefaultValuePluginFunction&
        ) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
