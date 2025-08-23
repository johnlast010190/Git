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
    Foam::findShiftedFieldPluginFunction

Description

SourceFiles
    findShiftedFieldPluginFunction.C

Contributors/Copyright:
    2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef findShiftedFieldPluginFunction_H
#define findShiftedFieldPluginFunction_H

#include "plugins/FieldValuePluginFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class findShiftedFieldPluginFunction Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class findShiftedFieldPluginFunction
:
    public FieldValuePluginFunction
{
protected:

    // Protected Data

        typedef VolField<Type> ResultType;

        Type defaultValue_;


    // Protected Member Functions

        void doEvaluation();

        void setArgument
        (
            label index,
            const string& content,
            const CommonValueExpressionDriver& driver
        );


private:

    // Private Data

        autoPtr<ResultType> field_;
        autoPtr<ResultType> default_;
        autoPtr<volVectorField> shift_;


public:

    //- Runtime type information
    TypeName("findShiftedFieldPluginFunction");


    // Constructors

        //- Construct from components
        findShiftedFieldPluginFunction
        (
            const FieldValueExpressionDriver& parentDriver,
            const word& name
        );

        //- Disallow default bitwise copy construction
        findShiftedFieldPluginFunction
        (
            const findShiftedFieldPluginFunction&
        ) = delete;


    //- Destructor
    virtual ~findShiftedFieldPluginFunction()
    {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const findShiftedFieldPluginFunction&) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
