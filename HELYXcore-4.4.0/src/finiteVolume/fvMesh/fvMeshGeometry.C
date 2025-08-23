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
    (c) 2011-2023 OpenFOAM Foundation
    (c) 2022-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMesh/fvMesh.H"
#include "db/Time/Time.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/volFields/slicedVolFields.H"
#include "fields/surfaceFields/slicedSurfaceFields.H"
#include "fields/Fields/Field/SubField.H"
#include "fields/fvPatchFields/constraint/cyclic/cyclicFvPatchFields.H"
#include "fields/fvPatchFields/constraint/cyclicAMI/cyclicAMIFvPatchFields.H"
#include "meshes/polyMesh/polyPatches/indirectPolyPatch/indirectPolyPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMesh::makeSf() const
{
    if (debug)
    {
        InfoInFunction << "Assembling face areas" << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (SfSlicePtr_ || SfPtr_)
    {
        FatalErrorInFunction
            << "face areas already exist"
            << abort(FatalError);
    }

    SfSlicePtr_ = new slicedSurfaceVectorField
    (
        IOobject
        (
            "Sf",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        *this,
        dimArea,
        faceAreas()
    );

    SfSlicePtr_->setOriented();

    slicedSurfaceVectorField::Boundary& SfB = SfSlicePtr_->boundaryFieldRef();

    forAll(boundary_, pI)
    {
        if (isA<indirectPolyPatch>(boundaryMesh()[pI]))
        {
            SfB[pI] = boundaryMesh()[pI].faceAreas();
        }
    }
}


void Foam::fvMesh::makeMagSf() const
{
    if (debug)
    {
        InfoInFunction << "Assembling mag face areas" << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (magSfSlicePtr_ || magSfPtr_)
    {
        FatalErrorInFunction
            << "mag face areas already exist"
            << abort(FatalError);
    }

    magSfSlicePtr_ = new slicedSurfaceScalarField
    (
        IOobject
        (
            "magSf",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        *this,
        dimArea,
        magFaceAreas()
    );
}


void Foam::fvMesh::makeC() const
{
    if (debug)
    {
        InfoInFunction << "Assembling cell centres" << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (CSlicePtr_ || CPtr_)
    {
        FatalErrorInFunction
            << "cell centres already exist"
            << abort(FatalError);
    }

    // Construct as slices. Only preserve processor (not e.g. cyclic)

    CSlicePtr_ = new slicedVolVectorField
    (
        IOobject
        (
            "Cc",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        *this,
        dimLength,
        cellCentres(),
        faceCentres(),
        true,               //preserveCouples
        true                //preserveProcOnly
    );
}


void Foam::fvMesh::makeCf() const
{
    if (debug)
    {
        InfoInFunction << "Assembling face centres" << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (CfSlicePtr_ || CfPtr_)
    {
        FatalErrorInFunction
            << "face centres already exist"
            << abort(FatalError);
    }

    tmp<vectorField> fc(faceCentres());

    CfSlicePtr_ = new slicedSurfaceVectorField
    (
        IOobject
        (
            "Cf",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        *this,
        dimLength,
        fc()
    );
}


Foam::surfaceVectorField& Foam::fvMesh::SfRef()
{
    if (!SfPtr_)
    {
        SfPtr_ = Sf().cloneUnSliced().ptr();

        deleteDemandDrivenData(SfSlicePtr_);
    }

    return *SfPtr_;
}


Foam::surfaceScalarField& Foam::fvMesh::magSfRef()
{
    if (!magSfPtr_)
    {
        magSfPtr_ = magSf().cloneUnSliced().ptr();

        deleteDemandDrivenData(magSfSlicePtr_);
    }

    return *magSfPtr_;
}


Foam::volVectorField& Foam::fvMesh::CRef()
{
    if (!CPtr_)
    {
        CPtr_ = C().cloneUnSliced().ptr();

        deleteDemandDrivenData(CSlicePtr_);
    }

    return *CPtr_;
}


Foam::surfaceVectorField& Foam::fvMesh::CfRef()
{
    if (!CfPtr_)
    {
        CfPtr_ = Cf().cloneUnSliced().ptr();

        deleteDemandDrivenData(CfSlicePtr_);
    }

    return *CfPtr_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volScalarField::Internal& Foam::fvMesh::V() const
{
    if (!VPtr_)
    {
        if (debug)
        {
            InfoInFunction
                << "Constructing from primitiveMesh::cellVolumes()" << endl;
        }

        VPtr_ = new SlicedDimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "Vc",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true
            ),
            *this,
            dimVolume,
            cellVolumes()
        );
    }

    return *VPtr_;
}


const Foam::volScalarField::Internal& Foam::fvMesh::V0() const
{
    if (!V0Ptr_ || Foam::isNull(V0Ptr_))
    {
        FatalErrorInFunction
            << "Vc0 is not available"
            << abort(FatalError);
    }

    return *V0Ptr_;
}


Foam::volScalarField::Internal& Foam::fvMesh::setV0()
{
    return const_cast<volScalarField::Internal&>(V0());
}


const Foam::volScalarField::Internal& Foam::fvMesh::V00() const
{
    if (!V00Ptr_)
    {
        V00Ptr_ = NullObjectPtr<DimensionedField<scalar, volMesh>>();
    }

    if (Foam::isNull(V00Ptr_))
    {
        return V0();
    }

    return *V00Ptr_;
}


Foam::tmp<Foam::volScalarField::Internal> Foam::fvMesh::Vsc() const
{
    if (moving() && time().subCycling())
    {
        const TimeState& ts = time();
        const TimeState& ts0 = time().prevTimeState();

        scalar tFrac =
        (
            ts.value() - (ts0.value() - ts0.deltaTValue())
        )/ts0.deltaTValue();

        if (tFrac < (1 - SMALL))
        {
            return V0() + tFrac*(V() - V0());
        }
        else
        {
            return V();
        }
    }
    else
    {
        return V();
    }
}


Foam::tmp<Foam::volScalarField::Internal> Foam::fvMesh::Vsc0() const
{
    if (moving() && time().subCycling())
    {
        const TimeState& ts = time();
        const TimeState& ts0 = time().prevTimeState();

        scalar t0Frac =
        (
            (ts.value() - ts.deltaTValue())
          - (ts0.value() - ts0.deltaTValue())
        )/ts0.deltaTValue();

        if (t0Frac > SMALL)
        {
            return V0() + t0Frac*(V() - V0());
        }
        else
        {
            return V0();
        }
    }
    else
    {
        return V0();
    }
}


const Foam::surfaceVectorField& Foam::fvMesh::Sf() const
{
    if (SfPtr_)
    {
        return *SfPtr_;
    }

    if (!SfSlicePtr_)
    {
        makeSf();
    }

    return *SfSlicePtr_;
}


const Foam::surfaceScalarField& Foam::fvMesh::magSf() const
{
    if (magSfPtr_)
    {
        return *magSfPtr_;
    }

    if (!magSfSlicePtr_)
    {
        makeMagSf();
    }

    return *magSfSlicePtr_;
}


const Foam::volVectorField& Foam::fvMesh::C() const
{
    if (CPtr_)
    {
        return *CPtr_;
    }

    if (!CSlicePtr_)
    {
        makeC();
    }

    return *CSlicePtr_;
}


const Foam::surfaceVectorField& Foam::fvMesh::Cf() const
{
    if (CfPtr_)
    {
        return *CfPtr_;
    }

    if (!CfSlicePtr_)
    {
        makeCf();
    }

    return *CfSlicePtr_;
}


Foam::tmp<Foam::surfaceVectorField> Foam::fvMesh::delta() const
{
    if (debug)
    {
        InfoInFunction << "Calculating face deltas" << endl;
    }

    tmp<surfaceVectorField> tdelta
    (
        surfaceVectorField::New("delta", *this, dimLength)
    );
    surfaceVectorField& delta = tdelta.ref();
    delta.setOriented();

    const volVectorField& C = this->C();
    const labelUList& owner = this->owner();
    const labelUList& neighbour = this->neighbour();

    forAll(owner, facei)
    {
        delta[facei] = C[neighbour[facei]] - C[owner[facei]];
    }

    surfaceVectorField::Boundary& deltabf =  delta.boundaryFieldRef();

    forAll(deltabf, patchi)
    {
        deltabf[patchi] = boundary()[patchi].delta();
    }

    return tdelta;
}


const Foam::surfaceScalarField& Foam::fvMesh::phi() const
{
    if (!phiPtr_)
    {
        FatalErrorInFunction
            << "mesh flux field does not exist, is the mesh actually moving?"
            << abort(FatalError);
    }

    // Set zero current time
    // mesh motion fluxes if the time has been incremented
    if (phiPtr_->timeIndex() != time().timeIndex())
    {
        (*phiPtr_) = dimensionedScalar(dimVolume/dimTime, 0);
    }

    phiPtr_->setOriented();

    return *phiPtr_;
}


Foam::surfaceScalarField& Foam::fvMesh::setPhi()
{
    if (!phiPtr_)
    {
        FatalErrorInFunction
            << "mesh flux field does not exist, is the mesh actually moving?"
            << abort(FatalError);
    }

    return *phiPtr_;
}


// ************************************************************************* //
