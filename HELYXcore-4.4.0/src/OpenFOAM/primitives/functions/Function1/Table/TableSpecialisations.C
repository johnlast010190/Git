/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
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
    (c) 2019-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/Table/TableBase.H"


// * * * * * * * * * * * * Template specialisations  * * * * * * * * * * * * //

template<>
bool Foam::Function1Types::TableBase<bool>::value
(
    const scalar x
) const
{
    scalar xDash = x;

    if (checkMinBounds(x, xDash))
    {
        return table_.first().second();
    }

    if (checkMaxBounds(xDash, xDash))
    {
        return table_.last().second();
    }

    // Use interpolator
    interpolator().valueWeights(xDash, currentIndices_, currentWeights_);

    for (label i = 0; i < currentIndices_.size(); i++)
    {
        if (currentWeights_[i] > 0)
        {
            return table_[currentIndices_[i]].second();
        }
    }

    WarningInFunction
        << "There are not positive weights: something is wrong"
        << endl;

    return false;
}


template<>
bool Foam::Function1Types::TableBase<bool>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    FatalErrorInFunction
        << "Table with bool type is not supported for integration"
        << exit(FatalError);
    return false;
}


template<>
bool Foam::Function1Types::TableBase<bool>::integrateYoverX
(
    const scalar x1,
    const scalar x2
) const
{
    FatalErrorInFunction
        << "Table with bool type is not supported for integration"
        << exit(FatalError);

    return false;
}


template<>
bool Foam::Function1Types::TableBase<bool>::derivative(const scalar x) const
{
    FatalErrorInFunction
        << "Table with bool type is not supported for derivative"
        << exit(FatalError);

    return false;
}


template<>
Foam::scalarField Foam::Function1Types::TableBase<Foam::scalarField>::derivative(const scalar x) const
{
    // Use interpolator
    interpolator().derivationWeights(x, currentIndices_, currentWeights_);
    scalarField sum(value(x).size(), pTraits<Foam::scalar>::zero);
    //Type sum = Zero;
    if (currentIndices_.size() == 1)
    {
        return sum;
    }
    for (label i = 0; i < currentIndices_.size() - 1; ++i)
    {
        sum +=
            (
                table_[currentIndices_[i + 1]].second()
                - table_[currentIndices_[i]].second()
            )/(currentWeights_[i + 1] - currentWeights_[i]);
    }
    if (currentIndices_.size() > 2)
    {
        sum /= (currentIndices_.size() - 1);
    }
    return sum;
}

// ************************************************************************* //
