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
    Foam::GeneralSetEvaluationPluginFunction

Description

SourceFiles
    GeneralSetEvaluationPluginFunction.C

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef GeneralSetEvaluationPluginFunction_H
#define GeneralSetEvaluationPluginFunction_H

#include "GeneralSetsPluginFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class GeneralSetEvaluationPluginFunction Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class GeneralSetEvaluationPluginFunction
:
    public GeneralSetsPluginFunction
{
    // Private Data

        autoPtr<Field<Type>> values_;


protected:

    // Protected Member Functions

        const Field<Type>& values() const
        {
            return values_();
        }

        void setArgument
        (
            label index,
            const string& content,
            const CommonValueExpressionDriver& driver
        );

        word getID(label index);

        typedef VolField<Type> resultType;


public:

    //- Runtime type information
    TypeName("GeneralSetEvaluationPluginFunction");


    // Constructors

        //- Construct from components
        GeneralSetEvaluationPluginFunction
        (
            const FieldValueExpressionDriver& parentDriver,
            const word& name
        );

        //- Disallow default bitwise copy construction
        GeneralSetEvaluationPluginFunction
        (
            const GeneralSetEvaluationPluginFunction&
        ) = delete;


    //- Destructor
    virtual ~GeneralSetEvaluationPluginFunction()
    {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const GeneralSetEvaluationPluginFunction&) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
