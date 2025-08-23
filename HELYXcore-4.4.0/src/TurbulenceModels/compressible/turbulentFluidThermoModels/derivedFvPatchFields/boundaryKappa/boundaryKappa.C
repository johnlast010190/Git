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
    (c) 2019-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "turbulentFluidThermoModels/derivedFvPatchFields/boundaryKappa/boundaryKappa.H"
#include "fields/volFields/volFields.H"
#include "fluidThermo/fluidThermo.H"
#include "solidThermo/solidThermo.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "multiphaseThermo/multiphaseThermo.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum<boundaryKappa::KMethodType, 3>::names[] =
        { "fluidThermo", "solidThermo", "lookup" };

    const NamedEnum<boundaryKappa::KMethodType, 3>
        boundaryKappa::KMethodTypeNames_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryKappa::boundaryKappa
(
    const objectRegistry& obr,
    const fvPatch& patch,
    const word& calculationType,
    const word& kappaName,
    const word& alphaAniName,
    const word& phaseName
)
:
    patch_(patch),
    obr_(obr),
    method_(KMethodTypeNames_[calculationType]),
    kappaName_(kappaName),
    alphaAniName_(alphaAniName),
    phaseName_(phaseName)
{}


Foam::boundaryKappa::boundaryKappa
(
    const objectRegistry& obr,
    const fvPatch& patch,
    const dictionary& dict,
    const word& phaseName
)
:
    patch_(patch),
    obr_(obr),
    method_(KMethodTypeNames_.read(dict.lookup("kappaMethod"))),
    kappaName_(dict.lookupOrDefault<word>("kappaName", "none")),
    alphaAniName_(dict.lookupOrDefault<word>("alphaAniName","Anialpha")),
    phaseName_(phaseName)
{}


Foam::boundaryKappa::boundaryKappa
(
    const objectRegistry& obr,
    const fvPatch& patch,
    const boundaryKappa& base
)
:
    patch_(patch),
    obr_(obr),
    method_(base.method_),
    kappaName_(base.kappaName_),
    alphaAniName_(base.alphaAniName_),
    phaseName_(base.phaseName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::boundaryKappa::nfAlphaEff() const
{
    const fvMesh& mesh = patch_.boundaryMesh().mesh();
    const label patchi = patch_.index();
    const word thermoName =
        IOobject::groupName(basicThermo::dictName, phaseName_);

    switch (method_)
    {
        case mtFluidThermo:
        {
            typedef compressible::turbulenceModel turbulenceModel;

            const word turbName =
                IOobject::groupName
                (
                    turbulenceModel::propertiesName,
                    phaseName_
                );

            if (obr_.foundObject<turbulenceModel>(turbName))
            {
                const turbulenceModel& turbModel =
                    obr_.lookupObject<turbulenceModel>(turbName);

                return patch_.nf()*turbModel.alphaEff(patchi);
            }
            else if (obr_.foundObject<basicThermo>(thermoName))
            {
                const basicThermo& thermo =
                    obr_.lookupObject<basicThermo>(thermoName);

                const scalarField& Tp = thermo.T().boundaryField()[patchi];
                scalarField CpByCpv
                (
                    thermo.Cp(Tp, patchi)/thermo.Cpv(Tp, patchi)
                );

                return
                    patch_.nf()
                   *(
                        thermo.kappa().boundaryField()[patchi]
                       /thermo.Cp().boundaryField()[patchi]
                    )*CpByCpv;
            }
            else
            {
                FatalErrorInFunction
                    << "Kappa defined to employ " << KMethodTypeNames_[method_]
                    << " method, but thermo package not available"
                    << exit(FatalError);
            }

            break;
        }

        case mtSolidThermo:
        {
            const basicThermo& bThermo =
                mesh.lookupObject<basicThermo>(thermoName);
            const solidThermo* solidThermoPtr = nullptr;
            if (isA<multiphaseThermo>(bThermo))
            {
                solidThermoPtr =
                    &dynamic_cast<const multiphaseThermo&>(bThermo).sThermo();
            }
            else
            {
                solidThermoPtr = dynamic_cast<const solidThermo*>(&bThermo);
            }
            const solidThermo& thermo = *solidThermoPtr;

            const scalarField& Tp = thermo.T().boundaryField()[patchi];
            scalarField CpByCpv
            (
                thermo.Cp(Tp, patchi)/thermo.Cpv(Tp, patchi)
            );

            if (thermo.isotropic())
            {
                return
                    patch_.nf()
                   *(
                        thermo.kappa().boundaryField()[patchi]
                       /thermo.Cp().boundaryField()[patchi]
                    )*CpByCpv;
            }
            else
            {
                const symmTensorField& alphaAni =
                    patch_.lookupPatchFieldInDb<volSymmTensorField, scalar>
                    (
                        obr_,
                        alphaAniName_
                    );

                return (patch_.nf() & alphaAni)*CpByCpv;
            }
            break;
        }

        case mtLookup:
        {
            const basicThermo& thermo =
                obr_.lookupObject<basicThermo>(thermoName);

            const scalarField& Tp = thermo.T().boundaryField()[patchi];
            scalarField Cpv(thermo.Cpv(Tp, patchi));
            return nfKappa()/Cpv;
        }

        default:
        {
            FatalErrorInFunction
                << "Unimplemented method " << KMethodTypeNames_[method_] << nl
                << "Please set 'kappaMethod' to one of "
                << flatOutput(KMethodTypeNames_.sortedToc()) << nl
                << "and 'kappaName' to the name of the volScalar"
                << " or volSymmTensor field (if kappaMethod=lookup)"
                << exit(FatalError);
        }
    }
    tmp<vectorField> tnull(new vectorField(0));

    return tnull;
}


Foam::tmp<Foam::vectorField> Foam::boundaryKappa::nfKappa() const
{
    const fvMesh& mesh = patch_.boundaryMesh().mesh();
    const label patchi = patch_.index();

    if (method_ == mtLookup)
    {
        if (obr_.foundObject<volScalarField>(kappaName_))
        {
            return
                patch_.nf()*patch_.lookupPatchFieldInDb<volScalarField, scalar>
                (
                    obr_,
                    kappaName_
                );
        }
        else if (obr_.foundObject<volSymmTensorField>(kappaName_))
        {
            const symmTensorField& KWall =
                patch_.lookupPatchFieldInDb<volSymmTensorField, scalar>
                (
                    obr_,
                    kappaName_
                );

            return (patch_.nf() & KWall);
        }
        else
        {
            FatalErrorInFunction
                << "Did not find field " << kappaName_
                << " on mesh " << mesh.name() << " patch " << patch_.name()
                << nl
                << "Please set 'kappaMethod' to one of "
                << flatOutput(KMethodTypeNames_.sortedToc()) << nl
                << "and 'kappa' to the name of the volScalar"
                << " or volSymmTensor field (if kappaMethod=lookup)"
                << exit(FatalError);

            tmp<vectorField> tnull(new vectorField(0));
            return tnull;
        }
    }
    else
    {
        const basicThermo& thermo =
            obr_.lookupObject<basicThermo>
            (
                IOobject::groupName
                (
                    basicThermo::dictName,
                    phaseName_
                )
            );
        const scalarField& Tp = thermo.T().boundaryField()[patchi];
        scalarField Cpv(thermo.Cp(Tp, patchi));
        return nfAlphaEff()*Cpv;
    }
}


Foam::tmp<Foam::scalarField> Foam::boundaryKappa::kappa() const
{
    return (nfKappa() & patch_.nf());
}


const Foam::basicThermo& Foam::boundaryKappa::thermo() const
{
    const word thermoName =
        IOobject::groupName(basicThermo::dictName, phaseName_);
    return obr_.lookupObject<basicThermo>(thermoName);
}


void Foam::boundaryKappa::write(Ostream& os) const
{
    os.writeEntry("kappaMethod", KMethodTypeNames_[method_]);
    if (method_ == mtLookup)
    {
        os.writeEntry("kappaName", kappaName_);
    }
    if (alphaAniName_ != "Anialpha")
    {
        os.writeEntry("alphaAniName", alphaAniName_);
    }
}


// ************************************************************************* //
