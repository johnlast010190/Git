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

Class
    Foam::expressionSource

Description
    Generates a field according to an expression

SourceFiles
    expressionSource.C

Contributors/Copyright:
    2010, 2012-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef expressionSource_H
#define expressionSource_H

#include "FieldValueExpressionDriver.H"
#include "include/swakFvOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class expressionSource Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class SwakBasicSourceCommon
:
    public swakFvOptionType
{
    // Private data

    SwakBasicSourceCommon(const SwakBasicSourceCommon&);

    FieldValueExpressionDriver driver_;

protected:

    List<exprString> expressions_;

    HashTable<dimensionSet,label> dimensions_;

    bool verbose_;

    FieldValueExpressionDriver &driver() { return driver_; }

    // Protected functions

        //- Set the local field data
        void setFieldData(const dictionary& dict, wordList& fieldNames);

public:

    //- Runtime type information
    //    TypeName("SwakBasicSourceCommon");

    // Constructors

        //- Construct from a dictionary
        SwakBasicSourceCommon
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

    // Destructor

        virtual ~SwakBasicSourceCommon();

    typedef VolField<T> resultField;

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames);

        // I-O

            //- Write the source properties
            virtual void writeData(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#   include "SwakBasicSourceCommon.C"
#endif

#endif

// ************************************************************************* //
