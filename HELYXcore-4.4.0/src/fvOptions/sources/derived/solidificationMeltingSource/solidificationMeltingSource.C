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
    (c) 2014-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "solidificationMeltingSource.H"
#include "fvMatrices/fvMatrices.H"
#include "basicThermo/basicThermo.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "fields/fvPatchFields/basic/extrapolatedCalculated/extrapolatedCalculatedFvPatchFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* NamedEnum
    <
        fv::solidificationMeltingSource::thermoMode,
        2
    >::names[] =
    {
        "thermo",
        "lookup"
    };

    namespace fv
    {
        defineTypeNameAndDebug(solidificationMeltingSource, 0);

        addToRunTimeSelectionTable
        (
            option,
            solidificationMeltingSource,
            dictionary
        );
    }
}

const Foam::NamedEnum<Foam::fv::solidificationMeltingSource::thermoMode, 2>
    Foam::fv::solidificationMeltingSource::thermoModeTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::fv::solidificationMeltingSource::Cp() const
{
    switch (mode_)
    {
        case mdThermo:
        {
            const basicThermo& thermo =
                obr_.lookupObject<basicThermo>(basicThermo::dictName);

            return volScalarField::New(thermo.Cp().name(), thermo.Cp());
            break;
        }
        case mdLookup:
        {
            if (CpName_ == "CpRef")
            {
                scalar CpRef = coeffs_.lookup<scalar>("CpRef");

                return tmp<volScalarField>
                (
                    volScalarField::New
                    (
                        name_ + "-Cp",
                        obr_,
                        mesh_,
                        dimensionedScalar
                        (
                            dimEnergy/dimMass/dimTemperature,
                            CpRef
                        ),
                        extrapolatedCalculatedFvPatchScalarField::typeName
                    )
                );
            }
            else
            {
                return obr_.lookupObject<volScalarField>(CpName_);
            }

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled thermo mode: " << thermoModeTypeNames_[mode_]
                << abort(FatalError);
        }
    }

    return tmp<volScalarField>(nullptr);
}


Foam::vector Foam::fv::solidificationMeltingSource::g() const
{
    if (obr_.foundObject<uniformDimensionedVectorField>("g"))
    {
        const uniformDimensionedVectorField& value =
            obr_.lookupObject<uniformDimensionedVectorField>("g");
        return value.value();
    }
    else
    {
        return coeffs_.lookup("g");
    }
}


void Foam::fv::solidificationMeltingSource::update(const volScalarField& Cp)
{
    if (curTimeIndex_ == mesh_.time().timeIndex())
    {
        return;
    }

    if (debug)
    {
        Info<< type() << ": " << name_ << " - updating phase indicator" << endl;
    }

    // update old time alpha1 field
    alpha1_.oldTime();

    const volScalarField& T = obr_.lookupObject<volScalarField>(TName_);

    forAll(cells_, i)
    {
        label celli = cells_[i];

        scalar Tc = T[celli];
        scalar Cpc = Cp[celli];
        scalar alpha1New = alpha1_[celli] + relax_*Cpc*(Tc - Tmelt_)/L_;

        alpha1_[celli] = max(0, min(alpha1New, 1));
        deltaT_[i] = Tc - Tmelt_;
    }

    alpha1_.correctBoundaryConditions();

    curTimeIndex_ = mesh_.time().timeIndex();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::solidificationMeltingSource::solidificationMeltingSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(sourceName, modelType, dict, obr),
    Tmelt_(coeffs_.lookup<scalar>("Tmelt")),
    L_(coeffs_.lookup<scalar>("L")),
    relax_(coeffs_.lookupOrDefault("relax", 0.9)),
    mode_(thermoModeTypeNames_.read(coeffs_.lookup("thermoMode"))),
    rhoRef_(coeffs_.lookup<scalar>("rhoRef")),
    TName_(coeffs_.lookupOrDefault<word>("T", "T")),
    CpName_(coeffs_.lookupOrDefault<word>("Cp", "Cp")),
    UName_(coeffs_.lookupOrDefault<word>("U", "U")),
    phiName_(coeffs_.lookupOrDefault<word>("phi", "phi")),
    Cu_(coeffs_.lookupOrDefault<scalar>("Cu", 100000)),
    q_(coeffs_.lookupOrDefault("q", 0.001)),
    beta_(coeffs_.lookup<scalar>("beta")),
    alpha1_
    (
        IOobject
        (
            name_ + ":alpha1",
            mesh_.time().timeName(),
            obr_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    curTimeIndex_(-1),
    deltaT_(cells_.size(), 0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::solidificationMeltingSource::sourceFields
(
    wordList& fieldNames
)
{
    fieldNames.setSize(2);
    fieldNames[0] = UName_;

    switch (mode_)
    {
        case mdThermo:
        {
            const basicThermo& thermo =
                obr_.lookupObject<basicThermo>(basicThermo::dictName);

            fieldNames[1] = thermo.T().name();
            break;
        }
        case mdLookup:
        {
            fieldNames[1] = TName_;
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled thermo mode: " << thermoModeTypeNames_[mode_]
                << abort(FatalError);
        }
    }
}


void Foam::fv::solidificationMeltingSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    apply(geometricOneField(), eqn);
}


void Foam::fv::solidificationMeltingSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    apply(rho, eqn);
}


void Foam::fv::solidificationMeltingSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    const volScalarField Cp(this->Cp());

    update(Cp);

    vector g = this->g();

    scalarField& Sp = eqn.diag();
    vectorField& Su = eqn.source();
    const scalarField& V = mesh_.V();

    forAll(cells_, i)
    {
        label celli = cells_[i];

        scalar Vc = V[celli];
        scalar alpha1c = alpha1_[celli];

        scalar S = -Cu_*sqr(1.0 - alpha1c)/(pow3(alpha1c) + q_);
        vector Sb = rhoRef_*g*beta_*deltaT_[i];

        Sp[celli] += Vc*S;
        Su[celli] += Vc*Sb;
    }
}


void Foam::fv::solidificationMeltingSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    // Momentum source uses a Boussinesq approximation - redirect
    addSup(eqn, fieldi);
}


// ************************************************************************* //
