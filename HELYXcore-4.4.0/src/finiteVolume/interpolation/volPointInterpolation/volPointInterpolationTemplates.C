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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "interpolation/volPointInterpolation/volPointInterpolation.H"
#include "interpolation/volPointInterpolation/pointConstraints.H"
#include "containers/Lists/CompactListList/CompactListList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::volPointInterpolation::interpolateUnconstrained
(
    const VolField<Type>& vf,
    PointField<Type>& pf
) const
{
    if (debug)
    {
        Pout<< "volPointInterpolation::interpolateUnconstrained("
            << "const VolField<Type>&, " << "PointField<Type>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

    const labelListList& pointCells = mesh().pointCells();
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    const fvBoundaryMesh& fvbm = mesh().boundary();

    const primitivePatch& boundary = boundaryPtr_();

    // Cache calls to patch coupled flags
    boolList isCoupledPolyPatch(pbm.size(), false);
    boolList isCoupledFvPatch(fvbm.size(), false);
    forAll(isCoupledFvPatch, patchi)
    {
        isCoupledPolyPatch[patchi] = pbm[patchi].coupled();
        isCoupledFvPatch[patchi] = fvbm[patchi].coupled();
    }

    pf.primitiveFieldRef() = Zero;
    pf.boundaryFieldRef() = Zero;

    // Interpolate from the cells
    forAll(pointWeights_, pointi)
    {
        forAll(pointWeights_[pointi], pointCelli)
        {
            const label celli = pointCells[pointi][pointCelli];

            pf[pointi] += pointWeights_[pointi][pointCelli]*vf[celli];
        }
    }

    // Get the boundary neighbour field
    const typename VolField<Type>::Boundary vfBnf
    (
        VolField<Type>::null(),
        vf.boundaryField().boundaryNeighbourField()
    );

    // Interpolate from the boundary faces
    forAll(boundaryPointWeights_, bPointi)
    {
        const label pointi = boundary.meshPoints()[bPointi];

        const labelList& pFaces = boundary.pointFaces()[bPointi];

        forAll(boundaryPointWeights_[bPointi], bPointFacei)
        {
            // FV indices
            const labelUList patches =
                mesh().polyBFacePatches()[pFaces[bPointFacei]];
            const labelUList patchFaces =
                mesh().polyBFacePatchFaces()[pFaces[bPointFacei]];

            forAll(boundaryPointWeights_[bPointi][bPointFacei], i)
            {
                // If FV coupled only, add the neighbouring cell value.
                if
                (
                    !isCoupledPolyPatch[patches[i]]
                 && isCoupledFvPatch[patches[i]]
                )
                {
                    pf[pointi] +=
                        boundaryPointNbrWeights_[bPointi][bPointFacei][i]
                       *vfBnf[patches[i]][patchFaces[i]];
                }

                // If not coupled, add a weight to the boundary value.
                if
                (
                    !isCoupledPolyPatch[patches[i]]
                 && !isCoupledFvPatch[patches[i]]
                 && !isA<indirectPolyPatch>(fvbm[patches[i]].patch())
                )
                {
                    pf[pointi] +=
                        boundaryPointWeights_[bPointi][bPointFacei][i]
                       *vf.boundaryField()[patches[i]][patchFaces[i]];
                }
            }
        }
    }

    // Synchronise
    syncTools::syncPointList(mesh(), pf, plusEqOp<Type>(), pTraits<Type>::zero);
}


template<class Type>
void Foam::volPointInterpolation::interpolate
(
    const VolField<Type>& vf,
    PointField<Type>& pf
) const
{
    interpolateUnconstrained(vf, pf);

    // Apply constraints
    pointConstraints::New(pf.mesh()).constrain(pf);
}


template<class Type>
Foam::tmp<Foam::PointField<Type>> Foam::volPointInterpolation::interpolate
(
    const VolField<Type>& vf,
    const word& name,
    const bool cache
) const
{
    typedef PointField<Type> PointFieldType;

    const pointMesh& pm = pointMesh::New(vf.mesh());
    const objectRegistry& db = pm.thisDb();

    if (!cache || vf.mesh().changing())
    {
        // Delete any old occurrences to avoid double registration
        if (db.objectRegistry::template foundObject<PointFieldType>(name))
        {
            PointFieldType& pf =
                db.objectRegistry::template lookupObjectRef<PointFieldType>
                (
                    name
                );

            if (pf.ownedByRegistry())
            {
                solution::cachePrintMessage("Deleting", name, vf);
                pf.release();
                delete &pf;
            }
        }

        tmp<PointField<Type>> tpf
        (
            PointField<Type>::New(name, pm, vf.dimensions())
        );

        interpolate(vf, tpf.ref());

        return tpf;
    }
    else
    {
        if (!db.objectRegistry::template foundObject<PointFieldType>(name))
        {
            solution::cachePrintMessage("Calculating and caching", name, vf);
            tmp<PointFieldType> tpf = interpolate(vf, name, false);
            PointFieldType* pfPtr = tpf.ptr();
            regIOobject::store(pfPtr);
            return *pfPtr;
        }
        else
        {
            PointFieldType& pf =
                db.objectRegistry::template lookupObjectRef<PointFieldType>
                (
                    name
                );

            if (pf.upToDate(vf))
            {
                solution::cachePrintMessage("Reusing", name, vf);
                return pf;
            }
            else
            {
                solution::cachePrintMessage("Deleting", name, vf);
                pf.release();
                delete &pf;

                solution::cachePrintMessage("Recalculating", name, vf);
                tmp<PointFieldType> tpf = interpolate(vf, name, false);

                solution::cachePrintMessage("Storing", name, vf);
                PointFieldType* pfPtr = tpf.ptr();
                regIOobject::store(pfPtr);

                return *pfPtr;
            }
        }
    }
}


template<class Type>
Foam::tmp<Foam::PointField<Type>>
Foam::volPointInterpolation::interpolate(const VolField<Type>& vf) const
{
    return interpolate(vf, "volPointInterpolate(" + vf.name() + ')', false);
}


template<class Type>
Foam::tmp<Foam::PointField<Type>> Foam::volPointInterpolation::interpolate
(
    const tmp<VolField<Type>>& tvf
) const
{
    tmp<PointField<Type>> tpf = interpolate(tvf());

    tvf.clear();

    return tpf;
}


template<class Type>
void Foam::volPointInterpolation::interpolateDimensionedInternalField
(
    const DimensionedField<Type, volMesh>& vf,
    DimensionedField<Type, pointMesh>& pf
) const
{
    if (debug)
    {
        Pout<< "volPointInterpolation::interpolateDimensionedInternalField("
            << "const DimensionedField<Type, volMesh>&, "
            << "DimensionedField<Type, pointMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

    const fvMesh& mesh = vf.mesh();

    const labelListList& pointCells = mesh.pointCells();
    const pointField& points = mesh.points();
    const vectorField& cellCentres = mesh.cellCentres();

    // Re-do weights and interpolation since normal interpolation
    // pointWeights_ are for non-boundary points only. Not efficient but
    // then saves on space.

    // Multiply volField by weighting factor matrix to create pointField
    scalarField sumW(points.size(), 0.0);
    forAll(pointCells, pointi)
    {
        const labelList& ppc = pointCells[pointi];

        pf[pointi] = Type(Zero);

        forAll(ppc, pointCelli)
        {
            label celli = ppc[pointCelli];
            scalar pw = 1.0/mag(points[pointi] - cellCentres[celli]);

            pf[pointi] += pw*vf[celli];
            sumW[pointi] += pw;
        }
    }

    // Sum collocated contributions
    pointConstraints::syncUntransformedData(mesh, sumW, plusEqOp<scalar>());
    pointConstraints::syncUntransformedData(mesh, pf, plusEqOp<Type>());

    // Normalise
    forAll(pf, pointi)
    {
        scalar s = sumW[pointi];
        if (s > ROOTVSMALL)
        {
            pf[pointi] /= s;
        }
    }
}


template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::pointMesh>>
Foam::volPointInterpolation::interpolate
(
    const DimensionedField<Type, volMesh>& vf
) const
{
    const word name("volPointInterpolate(" + vf.name() + ')');
    const bool cache = false;

    typedef DimensionedField<Type, pointMesh> PointFieldType;

    const pointMesh& pm = pointMesh::New(vf.mesh());
    const objectRegistry& db = pm.thisDb();

    if (!cache || vf.mesh().changing())
    {
        // Delete any old occurences to avoid double registration
        if (db.objectRegistry::template foundObject<PointFieldType>(name))
        {
            PointFieldType& pf = const_cast<PointFieldType&>
            (
                db.objectRegistry::template lookupObject<PointFieldType>(name)
            );

            if (pf.ownedByRegistry())
            {
                solution::cachePrintMessage("Deleting", name, vf);
                pf.release();
                delete &pf;
            }
        }

        tmp<DimensionedField<Type, pointMesh>> tpf
        (
            new DimensionedField<Type, pointMesh>
            (
                IOobject
                (
                    name,
                    vf.instance(),
                    pm.thisDb()
                ),
                pm,
                vf.dimensions()
            )
        );

        interpolateDimensionedInternalField(vf, tpf.ref());

        return tpf;
    }
    else
    {
        if (!db.objectRegistry::template foundObject<PointFieldType>(name))
        {
            solution::cachePrintMessage("Calculating and caching", name, vf);
            tmp<PointFieldType> tpf = interpolate(vf);
            PointFieldType* pfPtr = tpf.ptr();
            regIOobject::store(pfPtr);
            return *pfPtr;
        }
        else
        {
            PointFieldType& pf = const_cast<PointFieldType&>
            (
                db.objectRegistry::template lookupObject<PointFieldType>(name)
            );

            if (pf.upToDate(vf))    //TBD: , vf.mesh().points()))
            {
                solution::cachePrintMessage("Reusing", name, vf);
            }
            else
            {
                solution::cachePrintMessage("Updating", name, vf);
                interpolateDimensionedInternalField(vf, pf);
            }

            return pf;
        }
    }
}


// ************************************************************************* //
