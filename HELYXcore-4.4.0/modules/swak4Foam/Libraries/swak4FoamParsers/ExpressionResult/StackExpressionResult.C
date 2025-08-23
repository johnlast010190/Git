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
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "StackExpressionResult.H"
#include "primitives/Vector/vector/vector.H"
#include "primitives/Tensor/tensor/tensor.H"
#include "primitives/SymmTensor/symmTensor/symmTensor.H"
#include "primitives/SphericalTensor/sphericalTensor/sphericalTensor.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam {

defineTypeNameAndDebug(StackExpressionResult,0);

addToRunTimeSelectionTable(ExpressionResult, StackExpressionResult, dictionary);
addToRunTimeSelectionTable(ExpressionResult, StackExpressionResult, nothing);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

StackExpressionResult::StackExpressionResult()
:
    ExpressionResult()
{
    // this has to be reset every timestep to work
    setNeedsReset(true);
}

StackExpressionResult::StackExpressionResult(
    const StackExpressionResult &rhs
)
:
    ExpressionResult(rhs)
{
    // this has to be reset every timestep to work
    setNeedsReset(true);
}

StackExpressionResult::StackExpressionResult(const dictionary &dict)
:
    ExpressionResult(dict)
{
    // this has to be reset every timestep to work
    setNeedsReset(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

StackExpressionResult::~StackExpressionResult()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

ExpressionResult StackExpressionResult::pop()
{
    if (!(size()>0)) {
        FatalErrorIn("StackExpressionResult::pop()")
            << "Trying to pop result from a empty queue"
                << endl
                << abort(FatalError);

        return ExpressionResult();
    }

    if (valueType()==pTraits<scalar>::typeName) {
        return popInternal<scalar>();
    } else if (valueType()==pTraits<vector>::typeName) {
        return popInternal<vector>();
    } else if (valueType()==pTraits<tensor>::typeName) {
        return popInternal<tensor>();
    } else if (valueType()==pTraits<symmTensor>::typeName) {
        return popInternal<symmTensor>();
    } else if (valueType()==pTraits<sphericalTensor>::typeName) {
        return popInternal<sphericalTensor>();
    } else {
            FatalErrorIn("StackExpressionResult::pop()")
                << " Unsopported value type " << valueType()
                    << endl
                    << abort(FatalError);

            return ExpressionResult();
    }
}

void StackExpressionResult::push(ExpressionResult &atEnd)
{
    Dbug<< "push(ExpressionResult &atEnd)" << endl;
    Dbug<< "Pushing: " << atEnd << endl;

    if (!hasValue()) {
        // this is the first push
        //        static_cast<ExpressionResult>(*this)=atEnd;
        ExpressionResult::operator=(atEnd);
    } else {
        if (valueType()!=atEnd.valueType()) {
            FatalErrorIn("StackExpressionResult::push(const ExpressionResult &atEnd)")
                << "Type of pushed value " << atEnd.valueType()
                    << " is not the expected type " << valueType()
                    << endl
                    << abort(FatalError);
        }
        if (valueType()==pTraits<scalar>::typeName) {
            pushInternal<scalar>(atEnd);
        } else if (valueType()==pTraits<vector>::typeName) {
            pushInternal<vector>(atEnd);
        } else if (valueType()==pTraits<tensor>::typeName) {
            pushInternal<tensor>(atEnd);
        } else if (valueType()==pTraits<symmTensor>::typeName) {
            pushInternal<symmTensor>(atEnd);
        } else if (valueType()==pTraits<sphericalTensor>::typeName) {
            pushInternal<sphericalTensor>(atEnd);
        } else {
            FatalErrorIn("StackExpressionResult::push(const ExpressionResult &atEnd)")
                << " Unsopported value type " << valueType()
                    << endl
                    << abort(FatalError);
        }
    }

    Dbug<< "After push: " << *this << endl;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void StackExpressionResult::operator=(const StackExpressionResult& rhs)
{
    Dbug<< "operator=(const StackExpressionResult& rhs)"
        << endl;

    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("StackExpressionResult::operator=(const StackExpressionResult&)")
            << "Attempted assignment to self"
            << exit(FatalError);
    }

    static_cast<ExpressionResult&>(*this)=rhs;
}

void StackExpressionResult::operator=(const ExpressionResult& rhs)
{
    Dbug<< "operator=(const ExpressionResult& rhs)" << endl;

    ExpressionResult last(
        rhs.getUniform(
            1,
            false // issue a warning if the other result is not really uniform
        )
    );

    this->push(
        last
    );
}

} // namespace

// ************************************************************************* //
