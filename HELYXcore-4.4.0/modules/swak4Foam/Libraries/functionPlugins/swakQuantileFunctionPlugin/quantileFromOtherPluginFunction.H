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
    Foam::quantileFromOtherPluginFunction

Description

SourceFiles
    quantileFromOtherPluginFunction.C

Contributors/Copyright:
    2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef quantileFromOtherPluginFunction_H
#define quantileFromOtherPluginFunction_H

#include "plugins/FieldValuePluginFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class quantileFromOtherPluginFunction Declaration
\*---------------------------------------------------------------------------*/
class quantileFromOtherPluginFunction
:
    public FieldValuePluginFunction
{
    // Private Data

        autoPtr<volScalarField> frac_;
        autoPtr<volScalarField> other_;
        autoPtr<volScalarField> weight_;
        autoPtr<volScalarField> mask_;
        label binNumber_;


protected:

    // Protected Member Functions

        void setArgument
        (
            label index,
            const string& content,
            const CommonValueExpressionDriver& driver
        );

        void setArgument(label index, const label&);


public:

    //- Runtime type information
    TypeName("quantileFromOtherPluginFunction");


    // Constructors

        //- Construct from components
        quantileFromOtherPluginFunction
        (
            const FieldValueExpressionDriver& parentDriver,
            const word& name
        );

        //- Disallow default bitwise copy construction
        quantileFromOtherPluginFunction
        (
            const quantileFromOtherPluginFunction&
        ) = delete;


    //- Destructor
    virtual ~quantileFromOtherPluginFunction()
    {}


    // Member Functions

        void doEvaluation();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const quantileFromOtherPluginFunction&) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
