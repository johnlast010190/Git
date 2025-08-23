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
    (c) 2015 OpenCFD Ltd.
    (c) 2018-2020 OpenFOAM Foundation
    (c) 2021-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "discreteMixingPlanePolyPatch.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

bool Foam::discreteMixingPlanePolyPatch::is360() const
{
    return
    (
        (periodicPatchName_ == word::null)
      ||(periodicPatchName_ == "none")
    );
}

void Foam::discreteMixingPlanePolyPatch::syncTransforms() const
{
    if (owner())
    {
        // At this point we guarantee that the transformations have been
        // updated. There is one particular case, where if the periodic patch
        // locally has zero faces then its transformation will not be set. We
        // need to synchronise the transforms over the zero-sized patches as
        // well.
        //
        // We can't put the logic into cyclicPolyPatch as
        // processorCyclicPolyPatch uses cyclicPolyPatch functionality.
        // Synchronisation will fail because processor-type patches do not exist
        // on all processors.
        //
        // The use in cyclicPeriodicAMI is special; we use the patch even
        // though we have no faces in it. Ideally, in future, the transformation
        // logic will be abstracted out, and we won't need a periodic patch
        // here. Until then, we can just fix the behaviour of the zero-sized
        // coupled patches here

        // Get the periodic patch
        const coupledPolyPatch& periodicPatch
        (
            refCast<const coupledPolyPatch>
            (
                boundaryMesh()[periodicPatchID()]
            )
        );

        // If there are any zero-sized periodic patches
        if (returnReduce((size() && !periodicPatch.size()), orOp<bool>()))
        {
            // Note that zero-sized patches will have zero-sized fields for the
            // separation vector, forward and reverse transforms. These need
            // replacing with the transformations from other processors.

            // Parallel in this context refers to a parallel transformation,
            // rather than a rotational transformation.

            // Note that a cyclic with zero faces is considered parallel so
            // explicitly check for that.
            bool isParallel =
            (
                periodicPatch.size()
             && !periodicPatch.transform().transforms()
            );
            reduce(isParallel, orOp<bool>());

            if (isParallel)
            {
                // Sync a list of separation vectors
                List<vector> sep(Pstream::nProcs());
                sep[Pstream::myProcNo()] = periodicPatch.transform().t();
                Pstream::allGatherList(sep);

                // If locally we have zero faces pick the first one that has a
                // separation vector
                if (!periodicPatch.size())
                {
                    forAll(sep, procI)
                    {
                        if (sep[procI].size())
                        {
                            const_cast<vector&>
                            (
                                periodicPatch.transform().t()
                            ) = sep[procI];

                            break;
                        }
                    }
                }
            }
            else
            {
                // Sync a list of forward and reverse transforms
                List<tensor> forwardT(Pstream::nProcs());
                forwardT[Pstream::myProcNo()] = periodicPatch.transform().T();
                Pstream::allGatherList(forwardT);

                // If locally we have zero faces pick the first one that has a
                // transformation vector
                if (!periodicPatch.size())
                {
                    forAll(forwardT, procI)
                    {
                        if (forwardT[procI].size())
                        {
                            const_cast<tensor&>
                            (
                                periodicPatch.transform().T()
                            ) = forwardT[procI];

                            break;
                        }
                    }
                }
            }
        }
    }
}



void Foam::discreteMixingPlanePolyPatch::computeSector() const
{
    if (nSectors_ == 0)
    {
        if (sectorDefinition_ == periodicPatch)
        {
            if (periodicPatchID()==-1)
            {
                //- Corresponds to 360 | periodicPatchName_ == word::null
                nSectors_ = 1;
            }
            else
            {
                const cyclicTransform& cycTransfPatch
                (
                    dynamic_cast<const cyclicTransform&>
                    (
                        boundaryMesh()[periodicPatchID()]
                    )
                );

                const coupledPolyPatch& cpp
                (
                    refCast<const coupledPolyPatch>
                    (
                        boundaryMesh()[periodicPatchID()]
                    )
                );

                if (cycTransfPatch.transformType() != ROTATIONAL)
                {
                    FatalErrorInFunction
                        << "Periodic patch of cyclicPeriodicAMI patch "
                        << this->name() << " with name " << cpp.name()
                        << "has " << cycTransfPatch.transformType()
                        << " transformation type. " << endl
                        << "Currently only rotational is supported for "
                        << "cyclicPeriodicAMI patches."
                        << exit(FatalError);
                }

                scalar sec = cpp.sectors();
                scalar nSec = round(sec);

                if (mag(nSec-sec) > sectorMatchTolerance_)
                {
                    FatalErrorInFunction
                        << "Sector numbers are not computed correctly "
                        << "for patch name " << this->name()
                        << ", sector = "
                        << sec
                        << this->name()
                        << endl
                        << "Check mesh generation or if sector numbers is close"
                        << " to integer, then increase match tolerance: "
                        << matchTolerance()
                        << exit(FatalError);
                }

                nSectors_ = nSec;

            }

            if (nSectors_ == 1)
            {
                Info<< "Singe sector (360) for patch: "
                     << this->name()
                     << ", sectors = "
                     << nSectors_
                     << endl;
                if (nSubdivisions_ != 1)
                {
                    FatalErrorInFunction
                        << "Single sector can't use subdivisions for "
                        << "periodicPatch type definition "
                        << "for patch name " << this->name()
                        << " since the rotation can't be calculated." << nl
                        << " Use userDefined definition instead or "
                        << "change the nSubdivisions entry to 1"
                        << exit(FatalError);
                }
            }
            else
            {
                Info<< "Sectors not defined. "
                     << "Computed automatically from periodic patch "
                     << boundaryMesh()[periodicPatchID()].name()
                     << " for patch: " << this->name()
                     << ", sectors = " << nSectors_
                     << ", subdivisions = " << nSubdivisions_
                     << endl;
            }
        }
        else if (sectorDefinition_ == userDefined)
        {
            if (nSubdivisions_ == 0)
            {
                FatalErrorInFunction
                    << "nSectors or nSubdivisions not defined for patch." << nl
                    << "For userDefined sectorDefinition "
                    << "nSectors or nSubdivisions entry is required"
                     << this->name() << endl;
            }
        }
        else
        {
            FatalErrorInFunction
                << "sectorDefinition: userDefined or periodicPatch";
        }
    }
    else
    {
        if (sectorDefinition_ == userDefined)
        {
            Info<< "Averaging snapshots (nSectors) defined for patch "
                 << this->name()
                 << ", sectors = " << nSectors_
                 << ", subdivisions = " << nSubdivisions_
                 << endl;
        }
    }
}

Foam::label Foam::discreteMixingPlanePolyPatch::nSectors() const
{
    if (nSectors_==0)
    {
        computeSector();
    }
    if (nSectors_<1)
    {
        FatalErrorInFunction
            << "Sector for patch " << this->name()
            << " is not computed correctly"
            << endl
            << "nSectors: "
            << nSectors_
            << exit(FatalError);
    }
    const label totalSectors = nSectors_*nSubdivisions_;
    return totalSectors;
}


Foam::label Foam::discreteMixingPlanePolyPatch::thisSectors() const
{
    return this->nSectors();
}


Foam::label Foam::discreteMixingPlanePolyPatch::nbrSectors() const
{
    const discreteMixingPlanePolyPatch& cAMIpp =
        refCast<const discreteMixingPlanePolyPatch>
        (
            nbrPatch()
        );
    const label nbrSectors = cAMIpp.nSectors();

    return nbrSectors;
}



Foam::label Foam::discreteMixingPlanePolyPatch::periodicPatchID() const
{
    if (is360())
    {
        //- corresponds in 360 sector - no periodic patch name defined

        periodicPatchID_ = -1;

        return periodicPatchID_;
    }

    if (periodicPatchID_ == -1)
    {
        periodicPatchID_ = this->boundaryMesh().findPatchID(periodicPatchName_);

        if (periodicPatchID_ == -1)
        {
            FatalErrorInFunction
                << "Illegal periodicPatch name " << periodicPatchName_
                << nl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        // Check that it is a coupled patch
        refCast<const coupledPolyPatch>
        (
            this->boundaryMesh()[periodicPatchID_]
        );
    }

    return periodicPatchID_;
}


Foam::label Foam::discreteMixingPlanePolyPatch::nbrPeriodicPatchID() const
{
    const discreteMixingPlanePolyPatch& cAMIpp =
        refCast<const discreteMixingPlanePolyPatch>
        (
            nbrPatch()
        );

    return cAMIpp.periodicPatchID();
}


const Foam::transformer&
Foam::discreteMixingPlanePolyPatch::thisPatchTranformer() const
{
    if (thisTransformer_ == transformer::null)
    {
        if (sectorDefinition_ == periodicPatch)
        {
            const cyclicTransform& onCycTransform
            (
                refCast<const cyclicTransform>
                (
                    boundaryMesh()[periodicPatchID()]
                )
            );

            sectorRotationCentre_ = onCycTransform.rotationCentre();
            sectorRotationAxis_ = onCycTransform.rotationAxis();
        }
        const tensor R =
            quaternion
            (
                sectorRotationAxis_,
                2*constant::mathematical::pi/thisSectors()
            ).R();

        if (mag(sectorRotationCentre_) == 0)
        {
            thisTransformer_ = transformer::rotation(R);
        }
        else
        {
            thisTransformer_ =
            transformer::translation(sectorRotationCentre_)
            & transformer::rotation(R)
            & transformer::translation(-sectorRotationCentre_);
        }
    }
    return thisTransformer_;
}


const Foam::transformer&
Foam::discreteMixingPlanePolyPatch::nbrPatchTranformer() const
{
    if (nbrTransformer_ == transformer::null)
    {
        const discreteMixingPlanePolyPatch& cAMIpp =
            refCast<const discreteMixingPlanePolyPatch>
            (
                nbrPatch()
            );
        nbrTransformer_ = cAMIpp.thisPatchTranformer();
    }

    return nbrTransformer_;
}


Foam::transformer Foam::discreteMixingPlanePolyPatch::getTransform
(
    const coupledPolyPatch& transformPatch
) const
{
    // Get the transform associated with the transform patch
    transformer t;
    {
        Tuple2<bool, transformer> bt
        (
            transformPatch.size(),
            transformPatch.transform()
        );

        reduce(bt, keepIfTrueOp<transformer>());

        if (!bt.first())
        {
            FatalErrorInFunction
                << "Transform patch " << transformPatch.name() << " for "
                << typeName << " patch " << name() << " has zero faces. It may "
                << "have been converted to a processor cyclic during "
                << "decomposition. Consider adding " << transformPatch.name()
                << " and it's neighbour to the list of preserved patches."
                << exit(FatalError);
        }

        t = bt.second();
    }

    return t;
}


void Foam::discreteMixingPlanePolyPatch::writeOBJ
(
    const primitivePatch& p,
    OBJstream& str
) const
{
    // Collect faces and points
    pointField allPoints;
    faceList allFaces;
    labelList pointMergeMap;
    PatchTools::gatherAndMerge
    (
        -1.0,           // do not merge points
        p,
        allPoints,
        allFaces,
        pointMergeMap
    );

    if (Pstream::master())
    {
        // Write base geometry
        str.write(allFaces, allPoints);
    }
}


// ************************************************************************* //
