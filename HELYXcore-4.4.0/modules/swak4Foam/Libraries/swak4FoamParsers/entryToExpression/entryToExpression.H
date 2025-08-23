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

Class
    Foam::entryToExpression

Description
    Helper class that gets an entry and converts it into a valid
    swak-expression

SourceFiles
    entryToExpression.C

Contributors/Copyright:
    2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef entryToExpression_H
#define entryToExpression_H

#include "include/swak.H"

#include "db/dictionary/entry/entry.H"

#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class entryToExpression Declaration
\*---------------------------------------------------------------------------*/

class entryToExpression
{
    // Private data

    //- Construct as copy
    entryToExpression(const entryToExpression&);

    void operator=(const entryToExpression&);

protected:

    template<class Type>
    string toStr(
        const Type &data,
        const word &prefix
    );

public:
    //- Runtime type information
    TypeName("entryToExpression");

    declareRunTimeSelectionTable
    (
        autoPtr,
        entryToExpression,
        nothing,
        (
        ),
        ()
    );

    // Constructors

    //- Construct null
    entryToExpression();

   //- Return a reference to the selected value driver
    static autoPtr<entryToExpression> New
    (
        const word& name
    );

    // Destructor

    virtual ~entryToExpression();

    virtual string toExpr(const entry &) = 0;

    static string fromEntry(const entry &);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "entryToExpressionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
