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

Contributors/Copyright:
    2008-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "swakExpressionDistributionFunctionObject.H"
#include "fields/volFields/volFields.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <typename T>
void swakExpressionDistributionFunctionObject::getDistributionInternal(
    autoPtr<SimpleDistribution<T>> &dist,
    autoPtr<Field<T>> &sameWeight
) {
    if (sameWeight.valid()) {
        dist=setData(
            driver_->getResult<T>()(),
            sameWeight(),
            maskValues_()
        );
    } else if (weightValuesScalar_.valid()) {
        dist=setDataScalar(
            driver_->getResult<T>()(),
            weightValuesScalar_(),
            maskValues_()
        );
    } else {
        FatalErrorIn("swakExpressionDistributionFunctionObject::getDistributionInternal")
            << "Weight neither of type " << pTraits<scalar>::typeName
                << " nor " << pTraits<T>::typeName
                << endl
                << "Set weights are: "
                << pTraits<scalar>::typeName << ":"
                << weightValuesScalar_.valid() << " "
                << pTraits<vector>::typeName << ":"
                << weightValuesVector_.valid() << " "
                << pTraits<tensor>::typeName << ":"
                << weightValuesTensor_.valid() << " "
                << pTraits<symmTensor>::typeName << ":"
                << weightValuesSymmTensor_.valid() << " "
                << pTraits<sphericalTensor>::typeName << ":"
                << weightValuesSphericalTensor_.valid() << " "
                << exit(FatalError);

    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
