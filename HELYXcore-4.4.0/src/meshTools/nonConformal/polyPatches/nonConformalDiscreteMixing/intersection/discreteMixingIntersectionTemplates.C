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
    (c) 2024-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::discreteMixingIntersection::interpolateToNeighbour
(
    const Field<Type>& ownFld
) const
{
    tmp<Field<Type>> tFld(new Field<Type>(tgtAddress_.size(), Zero));
    Field<Type>& result = tFld.ref();

    // Define weighted sum operator
    const auto cop = multiplyWeightedOp<Type, plusEqOp<Type>>(plusEqOp<Type>());

    interpolateToNeighbour(ownFld, result, cop);

    return tFld;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::discreteMixingIntersection::interpolateToOwner
(
    const Field<Type>& nbrFld
) const
{
    tmp<Field<Type>> tFld(new Field<Type>(srcAddress_.size(), Zero));
    Field<Type>& result = tFld.ref();

    // Define weighted sum operator
    const auto cop = multiplyWeightedOp<Type, plusEqOp<Type>>(plusEqOp<Type>());

    interpolateToOwner(nbrFld, result, cop);

    return tFld;
}


template<class Type, class CombineOp>
void Foam::discreteMixingIntersection::interpolateToNeighbour
(
    const UList<Type>& ownFld,
    List<Type>& result,
    const CombineOp& cop
) const
{
    result.setSize(tgtAddress_.size());

    List<List<Type>> field(Pstream::nProcs());
    field[Pstream::myProcNo()] = ownFld;

    // Gather field if intersections are present across multiple processors
    if (Pstream::parRun())
    {
        Pstream::allGatherList(field);
    }

    const bool transformField = tgtIntersectionMap_.size();

    // Do the interpolation from owner to neighbour
    forAll(tgtAddress_, tgtFacei)
    {
        forAll(tgtAddress_[tgtFacei], i)
        {
            const remote& procFace = tgtAddress_[tgtFacei][i];
            Type value = field[procFace.proci][procFace.elementi];

            // Transform field, if necessary
            if (transformField)
            {
                const label ict = tgtIntersectionMap_[tgtFacei][i];
                const transformer& ownTransf = transforms_[ict].first();
                const transformer& nbrTransf = inv(transforms_[ict].second());

                value = ownTransf.transformVector(value);
                value = nbrTransf.transformVector(value);
            }

            // Add contributions to the interpolated field
            cop(result[tgtFacei], tgtFacei, value, tgtWeights_[tgtFacei][i]);
        }
    }

    // Stabilise the zero-weight faces. The value used for stabilisation
    // here is unimportant and just needs to be a sample of the interpolated
    // variable, since a zero area will be set for these faces anyway.
    if (tgtZeroWeightFaces_.size())
    {
        Type stabValue = result[0];
        forAll(tgtAddress_, tgtFacei)
        {
            if (!tgtZeroWeightFaces_.found(tgtFacei)) break;

            stabValue = result[tgtFacei + 1];
        }

        forAllConstIter(labelHashSet, tgtZeroWeightFaces_, iter)
        {
            cop(result[iter.key()], iter.key(), stabValue, 1.0);
        }
    }
}


template<class Type, class CombineOp>
void Foam::discreteMixingIntersection::interpolateToOwner
(
    const UList<Type>& nbrFld,
    List<Type>& result,
    const CombineOp& cop
) const
{
    result.setSize(srcAddress_.size());

    List<List<Type>> field(Pstream::nProcs());
    field[Pstream::myProcNo()] = nbrFld;

    // Gather field if intersections are present across multiple processors
    if (Pstream::parRun())
    {
        Pstream::allGatherList(field);
    }

    const bool transformField = srcIntersectionMap_.size();

    // Do the interpolation from neighbour to owner
    forAll(srcAddress_, srcFacei)
    {
        forAll(srcAddress_[srcFacei], i)
        {
            const remote& procFace = srcAddress_[srcFacei][i];
            Type value = field[procFace.proci][procFace.elementi];

            // Transform field, if necessary
            if (transformField)
            {
                const label ict = srcIntersectionMap_[srcFacei][i];
                const transformer& ownTransf = inv(transforms_[ict].first());
                const transformer& nbrTransf = transforms_[ict].second();

                value = nbrTransf.transformVector(value);
                value = ownTransf.transformVector(value);
            }

            // Add contributions to the interpolated field
            cop(result[srcFacei], srcFacei, value, srcWeights_[srcFacei][i]);
        }
    }

    // Stabilise the zero-weight faces. The value used for stabilisation
    // here is unimportant and just needs to be a sample of the interpolated
    // variable, since a zero area will be set for these faces anyway.
    if (srcZeroWeightFaces_.size())
    {
        Type stabValue = result[0];
        forAll(srcAddress_, srcFacei)
        {
            if (!srcZeroWeightFaces_.found(srcFacei)) break;

            stabValue = result[srcFacei + 1];
        }

        forAllConstIter(labelHashSet, srcZeroWeightFaces_, iter)
        {
            cop(result[iter.key()], iter.key(), stabValue, 1.0);
        }
    }
}


// ************************************************************************* //
