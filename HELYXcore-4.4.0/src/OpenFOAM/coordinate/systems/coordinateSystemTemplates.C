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
    (c) 2018 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "coordinate/systems/coordinateSystem.H"
#include "primitives/transform/transform.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class PointField>
Foam::tmp<Foam::tensorField>
Foam::coordinateSystem::rotationsImpl(const PointField& global) const
{
    const label len = global.size();

    auto tresult = tmp<tensorField>::New(len);
    auto& result = tresult.ref();

    for (label i=0; i<len; ++i)
    {
        result[i] = this->R(global[i]);
    }

    return tresult;
}


template<class PointField>
Foam::tmp<Foam::pointField>
Foam::coordinateSystem::transformPointImpl(const PointField& localCart) const
{
    const label len = localCart.size();

    auto tresult = tmp<pointField>::New(len);
    auto& result = tresult.ref();

    for (label i=0; i<len; ++i)
    {
        result[i] = Foam::transform(rot_, localCart[i]) + origin_;
    }

    return tresult;
}


template<class PointField>
Foam::tmp<Foam::pointField>
Foam::coordinateSystem::invTransformPointImpl(const PointField& global) const
{
    const label len = global.size();

    auto tresult = tmp<pointField>::New(len);
    auto& result = tresult.ref();

    for (label i=0; i<len; ++i)
    {
        result[i] = Foam::invTransform(rot_, global[i] - origin_);
    }

    return tresult;
}


template<class RetType, class Type, class BinaryOp>
Foam::tmp<Foam::Field<RetType>>
Foam::coordinateSystem::manyTimesImpl
(
    const tensor& tt,
    const UList<Type>& input,
    const BinaryOp& bop
)
{
    const label len = input.size();

    auto tresult = tmp<Field<RetType>>::New(len);
    auto& result = tresult.ref();

    for (label i=0; i<len; ++i)
    {
        result[i] = bop(tt, input[i]);
    }

    return tresult;
}


template<class RetType, class PointField, class Type, class BinaryOp>
Foam::tmp<Foam::Field<RetType>>
Foam::coordinateSystem::oneToOneImpl
(
    const PointField& global,
    const UList<Type>& input,
    const BinaryOp& bop
) const
{
    const label len = input.size();

    if (len != global.size())
    {
        FatalErrorInFunction
            << "positions has different size from input field"
            << abort(FatalError);
    }

    auto tresult = tmp<Field<RetType>>::New(len);
    auto& result = tresult.ref();

    for (label i=0; i<len; ++i)
    {
        result[i] = bop(this->R(global[i]), input[i]);
    }

    return tresult;
}


template<class RetType, class PointField, class Type, class BinaryOp>
Foam::tmp<Foam::Field<RetType>>
Foam::coordinateSystem::oneToManyImpl
(
    const PointField& global,
    const Type& input,
    const BinaryOp& bop
) const
{
    const label len = global.size();

    auto tresult = tmp<Field<RetType>>::New(len);
    auto& result = tresult.ref();

    for (label i=0; i<len; ++i)
    {
        result[i] = bop(this->R(global[i]), input);
    }

    return tresult;
}


// ************************************************************************* //
