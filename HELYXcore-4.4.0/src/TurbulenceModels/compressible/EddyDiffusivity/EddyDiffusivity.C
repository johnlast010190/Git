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
    (c) 2015-2019 OpenFOAM Foundation
    (c) 2020 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "EddyDiffusivity/EddyDiffusivity.H"
#include "turbulentFluidThermoModels/derivedFvPatchFields/wallFunctions/alphatWallFunctions/alphatKaderWallFunction/alphatKaderWallFunctionFvPatchScalarField.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void Foam::EddyDiffusivity<BasicTurbulenceModel>::correctNut()
{
    if (alphat_.valid())
    {
        alphat_() = this->rho_*this->nut()/Prt_;
        alphat_().correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::EddyDiffusivity<BasicTurbulenceModel>::EddyDiffusivity
(
    const word& type,
    const alphaField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicTurbulenceModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    // coeffsDict not available yet
    Prt_("Prt", dimless, 0.85),
    Sct_("Sct", dimless, 0.7),
    alphat_(nullptr)
{
    const word alphatName(IOobject::groupName("alphat", alphaRhoPhi.group()));
    tmp<volScalarField> alphaTemp(alphatReadIfPresent(alphatName));

    if (alphaTemp.valid())
    {
        alphat_.set(alphaTemp.ptr());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::EddyDiffusivity<BasicTurbulenceModel>::alphatReadIfPresent
(
    const word& fieldName
) const
{
    IOobject alphatHeader
    (
        fieldName,
        this->mesh_.time().timeName(),
        this->db_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (alphatHeader.typeHeaderOk<volScalarField>(true))
    {
        return tmp<volScalarField>
        (
            new volScalarField(alphatHeader, this->mesh_)
        );
    }
    else
    {
        return tmp<volScalarField>(nullptr);
    }
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::EddyDiffusivity<BasicTurbulenceModel>::alphatAutoCreate() const
{
    const word alphatName(IOobject::groupName("alphat", this->group()));
    Info<< "--> Creating " << alphatName
        << " to employ run-time selectable wall functions"
        << endl;

    const fvBoundaryMesh& bm = this->mesh_.boundary();

    wordList alphatBoundaryTypes(bm.size());

    forAll(bm, patchI)
    {
        if (isA<wallFvPatch>(bm[patchI]))
        {
            alphatBoundaryTypes[patchI] =
                compressible::
                alphatKaderWallFunctionFvPatchScalarField::typeName;
        }
        else
        {
            alphatBoundaryTypes[patchI] =
                calculatedFvPatchField<scalar>::typeName;
        }
    }

    tmp<volScalarField> alphat
    (
        new volScalarField
        (
            IOobject
            (
                alphatName,
                this->mesh_.time().timeName(),
                this->db_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            this->mesh_,
            dimensionedScalar(dimMass/(dimLength*dimTime), 0),
            alphatBoundaryTypes
        )
    );

    return alphat;
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::EddyDiffusivity<BasicTurbulenceModel>::alphat() const
{
    // Create alphat when needed
    if (!alphat_.valid())
    {
        alphat_.set(alphatAutoCreate().ptr());
        alphat_() = this->rho_*this->nut()/Prt_;
        alphat_().correctBoundaryConditions();
    }

    return alphat_();
}


template<class BasicTurbulenceModel>
bool Foam::EddyDiffusivity<BasicTurbulenceModel>::read()
{
    if (BasicTurbulenceModel::read())
    {
        Prt_.readIfPresent(this->coeffDict());
        Sct_.readIfPresent(this->coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::EddyDiffusivity<BasicTurbulenceModel>::Dmt() const
{
    tmp<volScalarField> Dmt = this->rho_*this->nut()/Sct_;
    Dmt->correctBoundaryConditions();

    return Dmt;
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::scalarField>
Foam::EddyDiffusivity<BasicTurbulenceModel>::Dmt(const label patchi) const
{
    tmp<scalarField> Dmt =
        this->rho_.boundaryField()[patchi]*this->nut(patchi)/Sct_.value();
    return Dmt;
}



template<class BasicTurbulenceModel>
void Foam::EddyDiffusivity<BasicTurbulenceModel>::correctEnergyTransport()
{
    EddyDiffusivity<BasicTurbulenceModel>::correctNut();
}


// ************************************************************************* //
