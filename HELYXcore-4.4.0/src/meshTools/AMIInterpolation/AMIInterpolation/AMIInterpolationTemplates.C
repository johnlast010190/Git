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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2015-2020 OpenCFD Ltd.
    (c) 2017-2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "AMIInterpolation/AMIInterpolation/AMIInterpolation.H"
#include "AMIInterpolation/AMIInterpolation/AMIMethod/AMIMethod/AMIMethod.H"
#include "meshTools/meshTools.H"
#include "meshes/polyMesh/polyDistributionMap/distributionMap.H"
#include "primitives/ops/flipOp.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#if defined(WIN64) || defined(WIN32)
#include "primitives/bools/Switch/Switch.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class CombineOp>
void Foam::AMIInterpolation::interpolateToTarget
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    if (fld.size() != srcAddress_.size())
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to source patch size" << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << "    supplied field = " << fld.size()
            << abort(FatalError);
    }

    if (lowWeightCorrection_ > 0)
    {
        if (defaultValues.size() != tgtAddress_.size())
        {
            FatalErrorInFunction
                << "Employing default values when sum of weights falls below "
                << lowWeightCorrection_
                << " but supplied default field size is not equal to target "
                << "patch size" << nl
                << "    default values = " << defaultValues.size() << nl
                << "    target patch   = " << tgtAddress_.size() << nl
                << abort(FatalError);
        }
    }

    result.setSize(tgtAddress_.size());

    bool transformed = (tgtAMITransforms_.size()>1);

    if (singlePatchProc_ == -1)
    {
        const distributionMap& map = srcMapPtr_();

        List<Type> work(fld);
        map.distribute(work);

        forAll(result, facei)
        {
            const labelList& faces = tgtAddress_[facei];
            const scalarList& weights = tgtWeights_[facei];
            scalar w2 = 1 - tgtWeightsSum_[facei];

            forAll(faces, i)
            {

                if (transformed)
                {
                    Type transformedType = work[faces[i]];

                    const labelList& transIDs = tgtTransforms_[facei];
                    const label transID = transIDs[i];
                    const transformer& tgtTrans = tgtAMITransforms_[transID];
                    transformedType = tgtTrans.transformVector(work[faces[i]]);

                    const transformer srcTrans =
                        inv(srcAMITransforms_[transID]);

                    transformedType = srcTrans.transformVector(transformedType);

                    cop(result[facei], facei, transformedType, weights[i]);
                }
                else
                {
                    cop(result[facei], facei, work[faces[i]], weights[i]);
                }
            }
            if (lowWeightCorrection_ > 0)
            {
                cop(result[facei], facei, defaultValues[facei], w2);
            }
        }
    }
    else
    {
        forAll(result, facei)
        {
            const labelList& faces = tgtAddress_[facei];
            const scalarList& weights = tgtWeights_[facei];
            scalar w2 = 1-tgtWeightsSum_[facei];

            forAll(faces, i)
            {
                if (transformed)
                {
                    Type transformedType = fld[faces[i]];

                    const labelList& transIDs = tgtTransforms_[facei];
                    const label transID = transIDs[i];
                    const transformer& tgtTrans = tgtAMITransforms_[transID];
                    transformedType = tgtTrans.transformVector(transformedType);

                    const transformer srcTrans =
                        inv(srcAMITransforms_[transID]);

                    transformedType = srcTrans.transformVector(transformedType);

                    cop(result[facei], facei, transformedType, weights[i]);
                }
                else
                {
                    cop(result[facei], facei, fld[faces[i]], weights[i]);
                }
            }
            if (lowWeightCorrection_ > 0)
            {
                cop(result[facei], facei, defaultValues[facei], w2);
            }
        }
    }
}


template<class Type, class CombineOp>
void Foam::AMIInterpolation::interpolateToSource
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    if (fld.size() != tgtAddress_.size())
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to target patch size" << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << "    supplied field = " << fld.size()
            << abort(FatalError);
    }

    if (lowWeightCorrection_ > 0)
    {
        if (defaultValues.size() != srcAddress_.size())
        {
            FatalErrorInFunction
                << "Employing default values when sum of weights falls below "
                << lowWeightCorrection_
                << " but supplied default field size is not equal to target "
                << "patch size" << nl
                << "    default values = " << defaultValues.size() << nl
                << "    source patch   = " << srcAddress_.size() << nl
                << abort(FatalError);
        }
    }

    result.setSize(srcAddress_.size());

    bool transformed = (srcAMITransforms_.size()>1);

    if (singlePatchProc_ == -1)
    {
        const distributionMap& map = tgtMapPtr_();

        List<Type> work(fld);
        map.distribute(work);


        forAll(result, facei)
        {
            const labelList& faces = srcAddress_[facei];
            const scalarList& weights = srcWeights_[facei];
            scalar w2 = 1-srcWeightsSum_[facei];
            forAll(faces, i)
            {
                Type transformedType = work[faces[i]];

                if (transformed)
                {
                    const labelList& transIDs = srcTransforms_[facei];
                    const label transID = transIDs[i];
                    const transformer& srcTrans = srcAMITransforms_[transID];
                    transformedType = srcTrans.transformVector(work[faces[i]]);

                    const transformer tgtTrans =
                        inv(tgtAMITransforms_[transID]);

                    transformedType = tgtTrans.transformVector(transformedType);

                    cop(result[facei], facei, transformedType, weights[i]);
                }
                else
                {
                    cop(result[facei], facei, work[faces[i]], weights[i]);
                }
            }
            if (lowWeightCorrection_ > 0)
            {
                cop(result[facei], facei, defaultValues[facei], w2);
            }
        }
    }
    else
    {
        forAll(result, facei)
        {
            const labelList& faces = srcAddress_[facei];
            const scalarList& weights = srcWeights_[facei];
            scalar w2 = 1-srcWeightsSum_[facei];

            forAll(faces, i)
            {
                if (transformed)
                {
                    Type transformedType = fld[faces[i]];

                    const labelList& transIDs = srcTransforms_[facei];
                    const label transID = transIDs[i];
                    const transformer& srcTrans = srcAMITransforms_[transID];
                    transformedType = srcTrans.transformVector(transformedType);

                    const transformer tgtTrans =
                        inv(tgtAMITransforms_[transID]);

                    transformedType = tgtTrans.transformVector(transformedType);

                    cop(result[facei], facei, transformedType, weights[i]);
                }
                else
                {
                    cop(result[facei], facei, fld[faces[i]], weights[i]);
                }
            }
            if (lowWeightCorrection_ > 0)
            {
                cop(result[facei], facei, defaultValues[facei], w2);
            }
        }
    }
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToSource
(
    const Field<Type>& fld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    tmp<Field<Type>> tresult
    (
        new Field<Type>
        (
            srcAddress_.size(),
            Zero
        )
    );

    interpolateToSource
    (
        fld,
        multiplyWeightedOp<Type, CombineOp>(cop),
        tresult.ref(),
        defaultValues
    );

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToSource
(
    const tmp<Field<Type>>& tFld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(tFld(), cop, defaultValues);
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToTarget
(
    const Field<Type>& fld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    tmp<Field<Type>> tresult
    (
        new Field<Type>
        (
            tgtAddress_.size(),
            Zero
        )
    );

    interpolateToTarget
    (
        fld,
        multiplyWeightedOp<Type, CombineOp>(cop),
        tresult.ref(),
        defaultValues
    );

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToTarget
(
    const tmp<Field<Type>>& tFld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(tFld(), cop, defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::AMIInterpolation::interpolateToSource
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(fld, plusEqOp<Type>(), defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToSource
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(tFld(), plusEqOp<Type>(), defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToTarget
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(fld, plusEqOp<Type>(), defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToTarget
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(tFld(), plusEqOp<Type>(), defaultValues);
}


// ************************************************************************* //
