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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2017-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMeshSubset/fvMeshSubset.H"
#include "fields/fvsPatchFields/constraint/internal/internalFvsPatchField.H"
#include "fields/pointPatchFields/constraint/internal/internalPointPatchField.H"
#include "fields/fvPatchFields/constraint/internal/internalFvPatchFields.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/directFvPatchFieldMapper.H"
#include "fvMeshAdder/directPointPatchFieldMapper.H"
#include "primitives/ops/flipOp.H"
#include "fvMesh/fvPatches/constraint/nonConformalCyclic/nonConformalCyclicFvPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::VolField<Type>>
Foam::fvMeshSubset::interpolate
(
    const VolField<Type>& vf,
    const fvMesh& sMesh,
    const labelList& patchMap,
    const labelList& cellMap,
    const labelList& faceMap,
    bool handleInternal
)
{
    // 1. Create the complete field with dummy patch fields
    PtrList<fvPatchField<Type>> patchFields(patchMap.size());

    forAll(patchFields, patchi)
    {
        // Set the first one by hand as it corresponds to the
        // exposed internal faces. Additional interpolation can be put here
        // as necessary.
        if (patchMap[patchi] == -1)
        {
            patchFields.set
            (
                patchi,
                new internalFvPatchField<Type>
                (
                    sMesh.boundary()[patchi],
                    DimensionedField<Type, volMesh>::null()
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchi,
                fvPatchField<Type>::New
                (
                    calculatedFvPatchField<Type>::typeName,
                    sMesh.boundary()[patchi],
                    DimensionedField<Type, volMesh>::null()
                )
            );
        }
    }

    tmp<VolField<Type>> tresF
    (
        new VolField<Type>
        (
            IOobject
            (
                "subset" + vf.name(),
                sMesh.time().timeName(),
                sMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sMesh,
            vf.dimensions(),
            Field<Type>(vf.primitiveField(), cellMap),
            patchFields
        )
    );
    VolField<Type>& resF = tresF.ref();
    resF.oriented() = vf.oriented();


    // 2. Change the fvPatchFields to the correct type using a mapper
    //  constructor (with reference to the now correct internal field)

    typename VolField<Type>::
        Boundary& bf = resF.boundaryFieldRef();

    forAll(bf, patchi)
    {
        if (patchMap[patchi] != -1)
        {
            // Construct addressing
            const fvPatch& basePatch = vf.mesh().boundary()[patchMap[patchi]];
            if (isA<nonConformalFvPatch>(basePatch)) continue;

            const label baseStart = basePatch.start();
            const label baseSize = basePatch.size();

            const fvPatch& subPatch = sMesh.boundary()[patchi];
            labelList directAddressing(subPatch.size());

            forAll(directAddressing, i)
            {
                const label baseFacei = faceMap[subPatch.start() + i];

                if (baseFacei >= baseStart && baseFacei < baseStart+baseSize)
                {
                    directAddressing[i] = baseFacei-baseStart;
                }
                else
                {
                    // Mapped from internal face. Do what? Leave up to
                    // fvPatchField
                    directAddressing[i] = -1;
                }
            }

            directFvPatchFieldMapper mapper(directAddressing);

            // suppress warning for unmapped faces (e.g. from fvMeshDistribute).
            if (handleInternal)
            {
                mapper.suppressWarning() = true;
            }

            bf.set
            (
                patchi,
                fvPatchField<Type>::New
                (
                    vf.boundaryField()[patchMap[patchi]],
                    subPatch,
                    resF(),
                    mapper
                )
            );
        }
    }

    if (handleInternal)
    {
        // Interpolate exposed internal faces.
        IStringStream schemeData("linear");
        SurfaceField<Type> vfs
        (
            fvc::interpolate(vf, schemeData)
        );

        forAll(bf, patchi)
        {
            if (patchMap[patchi] == -1)
            {
                fvPatchField<Type>& pfld = bf[patchi];

                label facei = pfld.patch().start();

                forAll(pfld, i)
                {
                    const label oldFacei = faceMap[facei++];
                    if (vf.mesh().isInternalFace(oldFacei))
                    {
                        pfld[i] = vfs[oldFacei];
                    }
                    else
                    {
                        const label mapPatchi =
                            vf.mesh().boundaryMesh().whichPatch(oldFacei);
                        const label index =
                            oldFacei - vf.mesh().boundary()[mapPatchi].start();

                        pfld[i] = vf.boundaryField()[mapPatchi][index];
                    }
                }
            }
        }
    }

    return tresF;
}


template<class Type>
Foam::tmp<Foam::VolField<Type>>
Foam::fvMeshSubset::interpolate
(
    const VolField<Type>& vf,
    bool handleInternal
) const
{
    return interpolate
    (
        vf,
        subMesh(),
        patchMap(),
        cellMap(),
        faceMap(),
        handleInternal
    );
}


template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::fvMeshSubset::interpolate
(
    const SurfaceField<Type>& sf,
    const fvMesh& sMesh,
    const labelList& patchMap,
    const labelList& cellMap,
    const labelList& faceMap,
    bool handleInternal
)
{
    // const bool negateIfFlipped = isFlux(sf);

    // 1. Create the complete field with dummy patch fields
    PtrList<fvsPatchField<Type>> patchFields(patchMap.size());

    forAll(patchFields, patchi)
    {
        // Set the first one by hand as it corresponds to the
        // exposed internal faces. Additional interpolation can be put here
        // as necessary.
        if (patchMap[patchi] == -1)
        {
            patchFields.set
            (
                patchi,
                new internalFvsPatchField<Type>
                (
                    sMesh.boundary()[patchi],
                    DimensionedField<Type, surfaceMesh>::null()
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchi,
                fvsPatchField<Type>::New
                (
                    calculatedFvsPatchField<Type>::typeName,
                    sMesh.boundary()[patchi],
                    DimensionedField<Type, surfaceMesh>::null()
                )
            );
        }
    }

    // Create the complete field from the pieces
    tmp<SurfaceField<Type>> tresF
    (
        new SurfaceField<Type>
        (
            IOobject
            (
                "subset" + sf.name(),
                sMesh.time().timeName(),
                sMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sMesh,
            sf.dimensions(),
            Field<Type>
            (
                sf.primitiveField(),
                SubList<label>
                (
                    faceMap,
                    sMesh.nInternalFaces()
                )
            ),
            patchFields
        )
    );
    SurfaceField<Type>& resF = tresF.ref();
    resF.oriented() = sf.oriented();


    // 2. Change the fvsPatchFields to the correct type using a mapper
    //  constructor (with reference to the now correct internal field)

    typename SurfaceField<Type>::
        Boundary& bf = resF.boundaryFieldRef();

    forAll(bf, patchi)
    {
        const fvPatch& subPatch = sMesh.boundary()[patchi];
        labelList directAddressing(subPatch.size(), -1);

        if (patchMap[patchi] != -1)
        {
            // Construct addressing
            const fvPatch& basePatch = sf.mesh().boundary()[patchMap[patchi]];
            const label baseStart = basePatch.start();
            const label baseSize = basePatch.size();

            forAll(directAddressing, i)
            {
                const label baseFacei = faceMap[subPatch.start() + i];

                if (baseFacei >= baseStart && baseFacei < baseStart+baseSize)
                {
                    directAddressing[i] = baseFacei-baseStart;
                }
            }

            bf.set
            (
                patchi,
                fvsPatchField<Type>::New
                (
                    sf.boundaryField()[patchMap[patchi]],
                    subPatch,
                    resF(),
                    directFvPatchFieldMapper(directAddressing)
                )
            );
        }

        // Map internal face values onto the patch elected to hold
        // the exposed faces
        fvsPatchField<Type>& pfld = bf[patchi];
        const labelUList& fc = bf[patchi].patch().faceCells();
        const labelList& own = sf.mesh().faceOwner();

        forAll(pfld, i)
        {
            const label baseFacei = faceMap[subPatch.start() + i];
            if (directAddressing[i] == -1)
            {
                if (sf.mesh().isInternalFace(baseFacei))
                {
                    const Type val = sf.internalField()[baseFacei];

                    if (cellMap[fc[i]] == own[baseFacei] || !sf.oriented()())
                    {
                        pfld[i] = val;
                    }
                    else
                    {
                        pfld[i] = flipOp()(val);
                    }
                }
                else
                {
                    const label basePatchi =
                        sf.mesh().boundaryMesh().patchID()
                        [baseFacei - sf.mesh().nInternalFaces()];
                    const label basePatchFacei =
                        sf.mesh().boundaryMesh()[basePatchi]
                       .whichFace(baseFacei);

                    pfld[i] = sf.boundaryField()[basePatchi][basePatchFacei];
                }
            }
        }
    }

    return tresF;
}


template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::fvMeshSubset::interpolate
(
    const SurfaceField<Type>& sf,
    bool handleInternal
) const
{
    return interpolate
    (
        sf,
        subMesh(),
        patchMap(),
        cellMap(),
        faceMap(),
        handleInternal
    );
}


template<class Type>
Foam::tmp<Foam::PointField<Type>>
Foam::fvMeshSubset::interpolate
(
    const PointField<Type>& vf,
    const pointMesh& sMesh,
    const labelList& patchMap,
    const labelList& pointMap
)
{
    if (debug)
    {
        Pout<< "fvMeshSubset::interpolate field " << vf.name() << endl;
    }
    // 1. Create the complete field with dummy patch fields
    PtrList<pointPatchField<Type>> patchFields(patchMap.size());

    forAll(patchFields, patchi)
    {
        // Set the first one by hand as it corresponds to the
        // exposed internal faces.  Additional interpolation can be put here
        // as necessary.
        if (patchMap[patchi] == -1)
        {
            patchFields.set
            (
                patchi,
                new internalPointPatchField<Type>
                (
                    sMesh.boundary()[patchi],
                    DimensionedField<Type, pointMesh>::null()
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchi,
                pointPatchField<Type>::New
                (
                    calculatedPointPatchField<Type>::typeName,
                    sMesh.boundary()[patchi],
                    DimensionedField<Type, pointMesh>::null()
                )
            );
        }
    }

    // Create the complete field from the pieces
    tmp<PointField<Type>> tresF
    (
        new PointField<Type>
        (
            IOobject
            (
                "subset" + vf.name(),
                sMesh.time().timeName(),
                sMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sMesh,
            vf.dimensions(),
            Field<Type>(vf.primitiveField(), pointMap),
            patchFields
        )
    );
    PointField<Type>& resF = tresF.ref();
    resF.oriented() = vf.oriented();


    // 2. Change the pointPatchFields to the correct type using a mapper
    //  constructor (with reference to the now correct internal field)

    typename PointField<Type>::
        Boundary& bf = resF.boundaryFieldRef();

    forAll(bf, patchi)
    {
        if (debug)
        {
            Pout<< "interpolate patch " << bf[patchi].patch().name()
                 << " for field " << resF.name() << endl;
        }
        // Set the first one by hand as it corresponds to the
        // exposed internal faces.  Additional interpolation can be put here
        // as necessary.
        if (patchMap[patchi] != -1)
        {
            // Construct addressing
            const pointPatch& basePatch =
                vf.mesh().boundary()[patchMap[patchi]];

            const labelList& meshPoints = basePatch.meshPoints();

            // Make addressing from mesh to patch point
            Map<label> meshPointMap(2*meshPoints.size());
            forAll(meshPoints, localI)
            {
                meshPointMap.insert(meshPoints[localI], localI);
            }

            // Find which subpatch points originate from which patch point
            const pointPatch& subPatch = sMesh.boundary()[patchi];
            const labelList& subMeshPoints = subPatch.meshPoints();

            // If mapped from outside patch leave handling up to patchField
            labelList directAddressing(subPatch.size(), -1);

            forAll(subMeshPoints, localI)
            {
                // Get mesh point on original mesh.
                label meshPointi = pointMap[subMeshPoints[localI]];

                Map<label>::const_iterator iter = meshPointMap.find(meshPointi);

                if (iter != meshPointMap.end())
                {
                    directAddressing[localI] = iter();
                }
            }

            bf.set
            (
                patchi,
                pointPatchField<Type>::New
                (
                    vf.boundaryField()[patchMap[patchi]],
                    subPatch,
                    resF(),
                    directPointPatchFieldMapper(directAddressing)
                )
            );
        }
    }

    return tresF;
}


template<class Type>
Foam::tmp<Foam::PointField<Type>>
Foam::fvMeshSubset::interpolate
(
    const PointField<Type>& sf,
    bool handleInternal
) const
{
    return interpolate
    (
        sf,
        pointMesh::New(subMesh()),     // subsetted point mesh
        patchMap(),
        pointMap()
    );
}


template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh>>
Foam::fvMeshSubset::interpolate
(
    const DimensionedField<Type, volMesh>& df,
    const fvMesh& sMesh,
    const labelList& cellMap
)
{
    // Create the complete field from the pieces
    tmp<DimensionedField<Type, volMesh>> tresF
    (
        new DimensionedField<Type, volMesh>
        (
            IOobject
            (
                "subset" + df.name(),
                sMesh.time().timeName(),
                sMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sMesh,
            df.dimensions(),
            Field<Type>(df, cellMap)
        )
    );

    tresF.ref().oriented() = df.oriented();

    return tresF;
}


template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh>>
Foam::fvMeshSubset::interpolate
(
    const DimensionedField<Type, volMesh>& df,
    bool handleInternal
) const
{
    return interpolate(df, subMesh(), cellMap());
}


// ************************************************************************* //
