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

Contributors/Copyright:
    2012-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2013 Georg Reiss <georg.reiss@ice-sf.at>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "StoredExpressionResult.H"
#include "primitives/Vector/vector/vector.H"
#include "primitives/Tensor/tensor/tensor.H"
#include "primitives/SymmTensor/symmTensor/symmTensor.H"
#include "primitives/SphericalTensor/sphericalTensor/sphericalTensor.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam {

defineTypeNameAndDebug(StoredExpressionResult,0);

addToRunTimeSelectionTable(ExpressionResult, StoredExpressionResult, dictionary);
addToRunTimeSelectionTable(ExpressionResult, StoredExpressionResult, nothing);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

StoredExpressionResult::StoredExpressionResult()
:
    ExpressionResult(),
    name_("None"),
    initialValueExpression_("")
{
}

StoredExpressionResult::StoredExpressionResult(
    const StoredExpressionResult &rhs
)
:
    ExpressionResult(rhs),
    name_(rhs.name_),
    initialValueExpression_(rhs.initialValueExpression_)
{
}

StoredExpressionResult::StoredExpressionResult(const dictionary &dict)
:
    ExpressionResult(dict.subOrEmptyDict("value")),
    name_(dict.lookup("name")),
    initialValueExpression_(
        dict.lookup("initialValue"),
        dict
    )
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

StoredExpressionResult::~StoredExpressionResult()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void StoredExpressionResult::operator=(const StoredExpressionResult& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("StoredExpressionResult::operator=(const StoredExpressionResult&)")
            << "Attempted assignment to self"
            << exit(FatalError);
    }

    static_cast<ExpressionResult&>(*this)=rhs;

    name_=rhs.name_;
    initialValueExpression_=rhs.initialValueExpression_;
}

void StoredExpressionResult::operator=(const ExpressionResult& rhs)
{
    //    static_cast<ExpressionResult&>(*this)=rhs;
   this->ExpressionResult::operator=(rhs);
}

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

Ostream & operator<<(Ostream &out,const StoredExpressionResult &data)
{
    out << token::BEGIN_BLOCK << endl;

    out.writeKeyword("name");
    out << word(data.name_) << token::END_STATEMENT << nl;

    out.writeKeyword("initialValue");
    out << data.initialValueExpression_ << token::END_STATEMENT << nl;

    out.writeKeyword("value");
    out << static_cast<const ExpressionResult &>(data);

    out << token::END_BLOCK << endl;

    return out;
}

Istream & operator>>(Istream &in,StoredExpressionResult &data)
{
    dictionary dict(in);

    data=StoredExpressionResult(dict);

    return in;
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

bool operator!=(const StoredExpressionResult &,const StoredExpressionResult &)
{
    return false;
}

} // namespace

// ************************************************************************* //
