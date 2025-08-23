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
    (c) 2020-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "addStaticHead/staticHead.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "faMesh/faMesh.H"
#include "include/faCFD.H"
#include "fields/faPatchFields/basic/zeroGradient/zeroGradientFaPatchField.H"
#include "fields/faPatchFields/basic/fixedGradient/fixedGradientFaPatchField.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"
#include "basicThermo/basicThermo.H"
#include "planarAverage/planarAverage.H"


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::staticHead::patchRho()
{
    basicThermo* thermoPtr =
        fvPf_.db().lookupObjectRefPtr<basicThermo>(basicThermo::matDictName);

    const fvPatch& fvP = fvPf_.patch();

    if (thermoPtr)
    {
        tmp<scalarField> prho
        (
            new scalarField
            (
                thermoPtr->distinctBuoyancy()
                ? thermoPtr->buoyantRho()().boundaryField()[fvP.index()]
                : thermoPtr->rho()().boundaryField()[fvP.index()]
            )
        );

        return prho;
    }

    thermoPtr =
        fvPf_.db().lookupObjectRefPtr<basicThermo>(basicThermo::dictName);

    if (thermoPtr)
    {
        return
            tmp<scalarField>
            (
                new scalarField(thermoPtr->rho()().boundaryField()[fvP.index()])
            );
    }

    // No thermo found - revert to rho lookup
    return
        tmp<scalarField>
        (
            new scalarField
            (
                fvPf_.lookupPatchField<volScalarField, scalar>(rhoName_)
            )
        );
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //


Foam::staticHead::staticHead
(
    const fvPatchField<scalar>& fvPf
)
:
    pStatic_(),
    fvPf_(fvPf),
    rhoRef_(1.0),
    useRhoAveAsRef_(false),
    rhoName_("rhoStaticHead"),
    integrateStaticHead_(false),
    hRefPoint_(nullptr),
    nNonOrthCorr_(0),
    addStaticHead_(false),
    nAveragingPlanes_(-1)
{}


Foam::staticHead::staticHead
(
    const fvPatchField<scalar>& fvPf,
    const fvPatch& p,
    const dictionary& dict
)
:
    pStatic_(),
    fvPf_(fvPf),
    rhoRef_(dict.lookupOrDefault<scalar>("rhoRef", 1.0)),
    useRhoAveAsRef_(!dict.found("rhoRef")),
    rhoName_(dict.lookupOrDefault<word>("rhoStaticHead", "rho")),
    integrateStaticHead_
    (
        dict.lookupOrDefault<Switch>("integrateStaticHead", false)
    ),
    hRefPoint_(nullptr),
    nNonOrthCorr_(dict.lookupOrDefault<label>("nNonOrthogonalCorrectors", 1)),
    addStaticHead_(dict.lookupOrDefault<Switch>("addStaticHead", false)),
    nAveragingPlanes_(dict.lookupOrDefault<label>("nAveragingPlanes", -1))
{
    if (integrateStaticHead_)
    {
        if (dict.found("pStatic"))
        {
            pStatic_ = scalarField("pStatic", dict, p.size());
        }
        else
        {
            pStatic_ = scalarField(p.size(), 0);
        }
    }
    if (dict.found("hRefPoint"))
    {
        hRefPoint_.set(new point(dict.lookup("hRefPoint")));
    }
}


Foam::staticHead::staticHead
(
    const fvPatchField<scalar>& fvPf,
    const staticHead& sH,
    const fvPatchFieldMapper& mapper
)
:
    pStatic_(mapper(sH.pStatic_)),
    fvPf_(fvPf),
    rhoRef_(sH.rhoRef_),
    useRhoAveAsRef_(sH.useRhoAveAsRef_),
    rhoName_(sH.rhoName_),
    integrateStaticHead_(sH.integrateStaticHead_),
    hRefPoint_(sH.hRefPoint_.valid() ? new point(sH.hRefPoint_()) : nullptr),
    nNonOrthCorr_(sH.nNonOrthCorr_),
    addStaticHead_(sH.addStaticHead_),
    nAveragingPlanes_(sH.nAveragingPlanes_)
{}


Foam::staticHead::staticHead
(
    const fvPatchField<scalar>& fvPf,
    const staticHead& sH
)
:
    pStatic_(sH.pStatic_),
    fvPf_(fvPf),
    rhoRef_(sH.rhoRef_),
    useRhoAveAsRef_(sH.useRhoAveAsRef_),
    rhoName_(sH.rhoName_),
    integrateStaticHead_(sH.integrateStaticHead_),
    hRefPoint_
    (
        sH.hRefPoint_.valid() ? new point(sH.hRefPoint_()) : nullptr
    ),
    nNonOrthCorr_(sH.nNonOrthCorr_),
    addStaticHead_(sH.addStaticHead_),
    nAveragingPlanes_(sH.nAveragingPlanes_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::staticHead::computeStaticHead()
{
    tmp<scalarField> tpsh(new scalarField(fvPf_.size(), 0));

    if (!active())
    {
        return tpsh;
    }

    scalarField& psh = tpsh.ref();

    dimensionedVector g("g", dimAcceleration, vector::zero);
    const uniformDimensionedVectorField* gPtr =
        fvPf_.db().lookupObjectPtr<uniformDimensionedVectorField>("g");

    const fvPatch& fvP = fvPf_.patch();

    if (gPtr)
    {
        g = *gPtr;
    }

    // Points adjacent to the cell centres, rather than face centres
    pointField patchFacePoints
    (
        fvP.patchInternalField(fvPf_.internalField().mesh().C())
      + fvP.delta()
    );

    scalarField gh
    (
        g.value() & patchFacePoints
    );

    scalar rhoRefMod = rhoRef_;
    if (useRhoAveAsRef_)
    {
        scalarField rhop(patchRho());
        rhoRefMod = gSum(rhop*fvP.magSf())/gSum(fvP.magSf());
    }

    // Read the solver hRef (if present) - used in p_rgh solvers
    dimensionedScalar hRef("hRef", dimLength, 0.0);
    const uniformDimensionedScalarField* hRefPtr =
        fvPf_.db().lookupObjectPtr<uniformDimensionedScalarField>("hRef");
    if (hRefPtr)
    {
        hRef = *hRefPtr;
    }

    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
        ? g & (cmptMag(g.value())/mag(g.value()))*hRef
        : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );

    if (integrateStaticHead_)
    {
        const fvMesh& mesh(fvPf_.internalField().mesh());

        // Stores the MeshObject for future reuse
        const faMesh& aMesh =
            MeshObject<polyMesh, UpdateableMeshObject, faMesh>::New
            (
                fvP.name(),
                mesh,
                fvP.name()
            );
        areaScalarField pPatch
        (
            IOobject
            (
                "p",
                mesh.time().timeName(),
                aMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh,
            dimensionedScalar(dimPressure, 0),
            fixedGradientFaPatchField<scalar>::typeName
        );
        pPatch.primitiveFieldRef() = pStatic_;

        scalarField rhop
        (
            planarAverage::planarAverage
            (
                aMesh,
                patchRho(),
                g.value()/stabilise(mag(g.value()), SMALL),
                nAveragingPlanes_
            )
        );

        areaScalarField rhoPatch
        (
            IOobject
            (
                "rho",
                mesh.time().timeName(),
                aMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh,
            dimensionedScalar(dimDensity, 0),
            zeroGradientFaPatchField<scalar>::typeName
        );
        rhoPatch.primitiveFieldRef() = rhop;
        rhoPatch.correctBoundaryConditions();

        forAll(pPatch.boundaryField(), bi)
        {
            if
            (
                isA<fixedGradientFaPatchField<scalar>>
                (
                    pPatch.boundaryField()[bi]
                )
            )
            {
                fixedGradientFaPatchField<scalar>& fgPatch =
                    refCast<fixedGradientFaPatchField<scalar>>
                    (
                        pPatch.boundaryFieldRef()[bi]
                    );
                fgPatch.gradient() =
                    (rhoPatch.boundaryField()[bi]*g.value())
                  & aMesh.boundary()[bi].edgeNormals();
            }
        }
        pPatch.correctBoundaryConditions();

        areaScalarField ghPatch(aMesh.areaCentres() & g);
        edgeScalarField ghfPatch(aMesh.edgeCentres() & g);

        word pName = fvPf_.internalField().name();

        // Ensure that we use the same snGrad scheme as the pressure laplacian,
        // for consistency
        tmp<fa::laplacianScheme<scalar, scalar>> pLaplacianScheme
        (
            fa::laplacianScheme<scalar, scalar>::New
            (
                aMesh,
                aMesh.schemes().laplacianScheme
                (
                    "laplacian("+pName+')'
                ),
                "grad("+pName+')'
            )
        );

        // Include p and subtract it again below to reduce clipping
        // by limited scheme (and be consistent with treatment in flow solver)
        edgeScalarField phig
        (
            "phig",
            (
                pLaplacianScheme().lnGradientScheme().lnGrad
                (
                    pPatch + rhoPatch*ghPatch
                )
              - pLaplacianScheme().lnGradientScheme().lnGrad(pPatch)
              - ghfPatch*pLaplacianScheme().lnGradientScheme().lnGrad(rhoPatch)
            )
           *aMesh.magLe()
        );

        // Find closest point if hRefPoint_ specified, or highest point if not
        scalarField diff;
        if (hRefPoint_.valid())
        {
            diff = mag(gh - (g.value() & hRefPoint_()));
        }
        else
        {
            diff = gh;
        }
        Tuple2<scalar,Pair<label>> minDiff
        (
            GREAT, Pair<label>(-1,-1)
        );
        forAll(diff, fi)
        {
            if (diff[fi] < minDiff.first())
            {
                minDiff.first() = diff[fi];
                minDiff.second().first() = fi;
            }
        }
        minDiff.second().second() = Pstream::myProcNo();
        reduce(minDiff, minAtOp<scalar>());

        if (hRefPoint_.valid())
        {
            // If hRefPoint_ was specified, its height must be within the
            // bounds of this patch
            scalar gMag = stabilise(mag(g.value()), SMALL);
            scalar minDist =
                minDiff.first()/gMag;
            scalar tol = 1e-6*mesh.bounds().mag();
            if (minDist > tol)
            {
                // Check with the patch points, not just face centres
                scalarField pgh
                (
                    g.value() & (fvP.patch().localPoints() - hRefPoint_())
                );
                scalar phMin = gMin(pgh)/gMag;
                scalar phMax = gMax(pgh)/gMag;
                if (phMin > tol || phMax < -tol)
                {
                    FatalErrorInFunction
                        << "Invalid height reference specified in patch "
                        << fvP.name()
                        << ":\nPlane containing hRefPoint_ with "
                        << "normal g must intersect some part of this patch, "
                        << "but its closest distance is " << max(phMin, -phMax)
                        << exit(FatalError);
                }
            }
        }

        label pRefFace(-1);
        scalar pRef(0);
        if (minDiff.second().second() == Pstream::myProcNo())
        {
            pRefFace = minDiff.second().first();
            // Set the pressure at the nearest face centre such as to be zero
            // at the actual point specified
            pRef =
                rhop[pRefFace]
               *(
                    (g.value() & fvPf_.patch().Cf()[pRefFace])
                  - (
                        hRefPoint_.valid()
                      ? (g.value() & hRefPoint_())
                        // Processor-local minimum is fine here as we know it'll
                        // be one of the points attached to pRefFace
                      : min(g.value() & fvP.patch().localPoints())
                    )
                );
        }

        for (label i = 0; i < nNonOrthCorr_+1; i++)
        {
            faMatrix<scalar> pEqn
            (
                -fam::laplacian(pPatch) + fac::div(phig)
            );
            pEqn.setReference(pRefFace, pRef);
            Info<< "Integrating static head on patch "
                << fvP.name() << ":"
                << endl;
            pEqn.solve();
        }
        pStatic_ = pPatch.primitiveField();
        // pStatic_ was calculated at the face centres, not at patchPoints
        scalarField pStaticCorr
        (
            pStatic_
          + (
                (patchFacePoints-fvP.Cf())
              & fac::grad(pPatch)().primitiveField()
            )
        );
        if (fvPf_.internalField().dimensions() == dimPressure/dimDensity)
        {
            psh += pStaticCorr/rhop;
        }
        else
        {
            psh += pStaticCorr;
        }
    }
    else
    {
        // The height reference used for calculating the static head
        scalar ghRefHead(0);
        if (hRefPoint_.valid())
        {
            // Use the specified value if given
            ghRefHead = (g.value() & hRefPoint_());
        }
        else
        {
            // Default to the solver hRef, or origin
            ghRefHead = ghRef.value();
        }
        psh += rhoRefMod*(gh - ghRefHead);
    }

    if (fvPf_.internalField().member() == "p_rgh")
    {
        scalarField rhop(patchRho());
        psh -= rhop*(gh - ghRef.value());
    }
    return tpsh;
}


void Foam::staticHead::addStaticHead(scalarField& pp)
{
    pp += computeStaticHead();
}


const Switch Foam::staticHead::active() const
{
    // Do not add static head if solving for potential flow, as gravity
    // is not accounted for
    if (fvPf_.internalField().name() == "pPot")
    {
            return false;
    }

    // With new thermo, only add static head if the buoyant() is true
    if (basicThermo::dictName == basicThermo::matDictName) // new thermo
    {
        basicThermo* thermoPtr =
            fvPf_.db().lookupObjectRefPtr<basicThermo>
            (
                basicThermo::matDictName
            );
        
        if (!thermoPtr) // no thermo
        {
            return addStaticHead_;
        }

        if (thermoPtr->buoyant())
        {
            return addStaticHead_;
        }
        else
        {
            return false;
        }
    }
    else // old thermo || no thermo
    {
        return addStaticHead_;
    }
}


void Foam::staticHead::autoMap
(
    const fvPatchFieldMapper& m
)
{
    if (integrateStaticHead_)
    {
        m(pStatic_, pStatic_);
    }
}


void Foam::staticHead::rmap
(
    const staticHead& sh, const labelList& addr
)
{
    if (integrateStaticHead_)
    {
        pStatic_.rmap(sh.pStatic_, addr);
    }
}


void Foam::staticHead::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    if (integrateStaticHead_)
    {
        mapper.map(pStatic_, scalar(0));
    }
}


void Foam::staticHead::write
(
    Ostream& os
)
const
{
    os.writeEntry("addStaticHead", addStaticHead_);
    if (addStaticHead_)
    {
        os.writeEntryIfDifferent("rhoStaticHead", word("rho"), rhoName_);
        if (!useRhoAveAsRef_)
        {
            os.writeEntry("rhoRef", rhoRef_);
        }
        os.writeEntry("integrateStaticHead", integrateStaticHead_);

        if (integrateStaticHead_)
        {
            pStatic_.writeEntry("pStatic", os);
            os.writeEntry("nNonOrthogonalCorrectors", nNonOrthCorr_);
            os.writeEntryIfDifferent
            (
                "nAveragingPlanes", label(-1), nAveragingPlanes_
            );
        }

        if (hRefPoint_.valid())
        {
            os.writeEntry("hRefPoint", hRefPoint_());
        }
    }
}


word Foam::staticHead::getRhoName() const
{
    return rhoName_;
}

// ************************************************************************* //
