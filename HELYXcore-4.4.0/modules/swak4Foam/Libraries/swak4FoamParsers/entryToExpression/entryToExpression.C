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

Contributors/Copyright:
    2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "entryToExpression.H"

#include "db/dictionary/primitiveEntry/primitiveEntry.H"
#include "db/IOstreams/StringStreams/OStringStream.H"

namespace Foam {

defineTypeNameAndDebug(entryToExpression,0);

defineRunTimeSelectionTable(entryToExpression, nothing);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

entryToExpression::entryToExpression()
{
}


autoPtr<entryToExpression> entryToExpression::New
(
    const word& name
)
{
    const auto ctor =
        ctorTableLookup
        (
            "entryToExpression type",
            nothingConstructorTable_(),
            name
        );
    return autoPtr<entryToExpression>(ctor());
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

entryToExpression::~entryToExpression()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

string entryToExpression::fromEntry(const entry &e)
{
    const primitiveEntry &pe=dynamicCast<const primitiveEntry&>(e);
    OStringStream o;

    for (label i=0; i<pe.size(); i++) {
        o << pe[i];

        if (i < pe.size()-1){
            o << token::SPACE;
        }
    }
    return o.str();
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace

// ************************************************************************* //
