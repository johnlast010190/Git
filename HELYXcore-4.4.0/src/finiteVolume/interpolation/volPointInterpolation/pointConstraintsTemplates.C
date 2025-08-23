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
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "interpolation/volPointInterpolation/pointConstraints.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "fields/pointPatchFields/basic/value/valuePointPatchFields.H"
#include "interpolation/volPointInterpolation/pointProcConstraintData/pointProcConstraintData.H"
#include "fields/pointPatchFields/constraint/processor/processorPointPatchField.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "primitives/transform/transformer/transformer.H"
#include "interpolation/volPointInterpolation/pointProcConstraintData/pointProcConstraintData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class CombineOp>
void Foam::pointConstraints::syncUntransformedData
(
    const polyMesh& mesh,
    List<Type>& pointData,
    const CombineOp& cop
)
{
    // Transfer onto coupled patch
    const globalMeshData& gmd = mesh.globalData();
    const indirectPrimitivePatch& cpp = gmd.coupledPatch();
    const labelList& meshPoints = cpp.meshPoints();

    const distributionMap& slavesMap = gmd.globalCoPointSlavesMap();
    const labelListList& slaves = gmd.globalCoPointSlaves();

    List<Type> elems(slavesMap.constructSize());
    forAll(meshPoints, i)
    {
        elems[i] = pointData[meshPoints[i]];
    }

    // Pull slave data onto master. No need to update transformed slots.
    slavesMap.distribute(elems, false);

    // Combine master data with slave data
    forAll(slaves, i)
    {
        Type& elem = elems[i];

        const labelList& slavePoints = slaves[i];

        // Combine master with untransformed slave data
        forAll(slavePoints, j)
        {
            cop(elem, elems[slavePoints[j]]);
        }

        // Copy result back to slave slots
        forAll(slavePoints, j)
        {
            elems[slavePoints[j]] = elem;
        }
    }

    // Push slave-slot data back to slaves
    slavesMap.reverseDistribute(elems.size(), elems, false);

    // Extract back onto mesh
    forAll(meshPoints, i)
    {
        pointData[meshPoints[i]] = elems[i];
    }
}


template<class Type>
void Foam::pointConstraints::setPatchFields
(
    PointField<Type>& pf
)
{
    typename PointField<Type>::
        Boundary& pfbf = pf.boundaryFieldRef();

    forAll(pfbf, patchi)
    {
        pointPatchField<Type>& ppf = pfbf[patchi];

        if (isA<valuePointPatchField<Type>>(ppf))
        {
            refCast<valuePointPatchField<Type>>(ppf) =
                ppf.patchInternalField();
        }
    }
}


template<class Type>
void Foam::pointConstraints::constrainCorners
(
    PointField<Type>& pf
) const
{
    forAll(patchPatchPointConstraintPoints_, pointi)
    {
        pf[patchPatchPointConstraintPoints_[pointi]] = transform
        (
            patchPatchPointConstraintTensors_[pointi],
            pf[patchPatchPointConstraintPoints_[pointi]]
        );
    }
}


template<class Type>
void Foam::pointConstraints::constrain
(
    PointField<Type>& pf,
    const bool overrideFixedValue
) const
{
    // Override constrained pointPatchField types with the constraint value.
    // This relies on only constrained pointPatchField implementing the evaluate
    // function
    pf.correctBoundaryConditions();

    if (!overrideFixedValue && Pstream::parRun())
    {
        constrainProcPoints(pf);
    }

    // Sync any dangling points
    syncUntransformedData
    (
        mesh()(),
        pf.primitiveFieldRef(),
        maxMagSqrEqOp<Type>()
    );

    // Apply multiple constraints on edge/corner points
    constrainCorners(pf);

    if (overrideFixedValue)
    {
        setPatchFields(pf);
    }
}


template<class Type>
void Foam::pointConstraints::constrainProcPoints
(
    PointField<Type>& pf
) const
{
    const pointField& pts = mesh()().points();
    List<pointProcConstraintData<Type>> lpData(pts.size());

    typename PointField<Type>::
        Boundary& pfbf = pf.boundaryFieldRef();

    //- Fill lpData
    //- mark proc points
    forAll(pfbf, pI)
    {
        pointPatchField<Type>& ppf = pfbf[pI];
        if (isA<processorPointPatchField<Type>>(ppf))
        {
            const labelList& ppp = ppf.patch().meshPoints();
            forAll(ppp, ppI)
            {
                const label gpI = ppp[ppI];
                lpData[gpI].procPoint() = true;
            }
        }
    }
    //- mark proc points that belong to other non-coupled boundaries
    forAll(pfbf, pI)
    {
        pointPatchField<Type>& ppf = pfbf[pI];
        if (!isA<coupledPointPatchField<Type>>(ppf))
        {
            const labelList& ppp = ppf.patch().meshPoints();
            forAll(ppp, ppI)
            {
                const label gpI = ppp[ppI];
                if (lpData[gpI].procPoint())
                {
                    lpData[gpI].procBoundaryPoint() = true;
                }
            }
        }
    }

    //- Prioritize non-fixed value. Take the corrected internal value
    forAll(pfbf, pI)
    {
        pointPatchField<Type>& ppf = pfbf[pI];

        if (isA<valuePointPatchField<Type>>(ppf))
        {
            // const valuePointPatchField<Type>& vppf =
            //    refCast<const valuePointPatchField<Type>>(ppf);
            Field<Type> valueField(ppf.patchInternalField()());
            const labelList& ppp = ppf.patch().meshPoints();
            forAll(ppp, ppI)
            {
                const label gpI = ppp[ppI];
                if (lpData[gpI].procBoundaryPoint())
                {
                    lpData[gpI].fixesValue() = true;
                    lpData[gpI].fixedValue() = valueField[ppI];
                }
            }
        }
    }

    //- sync
    syncUntransformedData
    (
        mesh()(),
        lpData,
        maxEqOp<pointProcConstraintData<Type>>()
    );

    //- fix the internal value
    forAll(lpData, gpI)
    {
        if (lpData[gpI].fixesValue())
        {
            pf.ref()[gpI] = lpData[gpI].fixedValue();
        }
    }
}


// ************************************************************************* //
