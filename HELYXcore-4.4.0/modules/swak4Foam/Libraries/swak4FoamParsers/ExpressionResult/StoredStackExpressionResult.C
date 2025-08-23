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

#include "StoredStackExpressionResult.H"
#include "primitives/Vector/vector/vector.H"
#include "primitives/Tensor/tensor/tensor.H"
#include "primitives/SymmTensor/symmTensor/symmTensor.H"
#include "primitives/SphericalTensor/sphericalTensor/sphericalTensor.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam {

defineTypeNameAndDebug(StoredStackExpressionResult,0);

addToRunTimeSelectionTable(ExpressionResult, StoredStackExpressionResult, dictionary);
addToRunTimeSelectionTable(ExpressionResult, StoredStackExpressionResult, nothing);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

StoredStackExpressionResult::StoredStackExpressionResult()
:
    StackExpressionResult()
{
    // reset the setting of the parent-class
    setNeedsReset(false);
}

StoredStackExpressionResult::StoredStackExpressionResult(
    const StoredStackExpressionResult &rhs
)
:
    StackExpressionResult(rhs)
{
    // reset the setting of the parent-class
    setNeedsReset(false);
}

StoredStackExpressionResult::StoredStackExpressionResult(const dictionary &dict)
:
    StackExpressionResult(dict)
{
    Dbug<< "StoredStackExpressionResult(const dictionary &dict)" << endl;
    Dbug<< "Value: " << (*this) << endl;

    // reset the setting of the parent-class
    setNeedsReset(false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

StoredStackExpressionResult::~StoredStackExpressionResult()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void StoredStackExpressionResult::resetInternal()
{
    // do nothing
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void StoredStackExpressionResult::operator=(const StoredStackExpressionResult& rhs)
{
    Dbug<< "operator=(const StoredStackExpressionResult& rhs)" << endl;

    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("StoredStackExpressionResult::operator=(const StoredStackExpressionResult&)")
            << "Attempted assignment to self"
            << exit(FatalError);
    }

    static_cast<StackExpressionResult&>(*this)=rhs;
}

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


} // namespace

// ************************************************************************* //
