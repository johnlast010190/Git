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
    (c) 2018-2024 Engys Ltd.
    (c) 2011-2022 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "absorptionEmissionSource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "radiationModels/fvDOM/fvDOM/fvDOM.H"
#include "global/constants/physicoChemical/physicoChemicalConstants.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "radiationModels/fvDOM/radiativeIntensityRay/radiativeIntensityRay.H"
#include "transportModel/transportModel.H"
#include "fluidThermo/fluidThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(absorptionEmissionSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        absorptionEmissionSource,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::absorptionEmissionSource::absorptionEmissionSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    TName_(coeffs_.lookupOrDefault<word>("TName", "T")),
    solarHeat_(coeffs_.lookupOrDefault<bool>("solarHeat", false)),
    e_(coeffs_.lookup<scalar>("emissivity")),
    a_(e_),
    E_(Function1<scalar>::New("E", coeffs_))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::absorptionEmissionSource::sourceFields(wordList& fieldNames)
{
    fieldNames = coeffs_.lookup<wordList>("fields");
}


void Foam::fv::absorptionEmissionSource::calculateRu(scalarField& Ru)
{
    const volScalarField& G_ = obr_.lookupObject<volScalarField>("G");

    forAll(cells_, i)
    {
        Ru[i] = a_*G_[cells_[i]] - E_->value(obr_.time().timeOutputValue());
    }

    if (solarHeat_ && obr_.foundObject<volScalarField>("Isolar"))
    {
        // absorbed fraction of the solar load dumped into heat source
        const volScalarField& Isolar_  = obr_.lookupObject<volScalarField>("Isolar");
        forAll(cells_, i)
        {
            Ru[i] += a_*Isolar_[cells_[i]];  // [W/m3]
        }
    }
}


void Foam::fv::absorptionEmissionSource::calculateRp(scalarField& Rp)
{
    forAll(cells_, i)
    {
        Rp[i] = 4.0*a_*constant::physicoChemical::sigma.value();
    }
}

void Foam::fv::absorptionEmissionSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{

    if (eqn.psi().name() == word(TName_))
    {
        // add source term on temperature equation
        addSupT(eqn);
    }

}


void Foam::fv::absorptionEmissionSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{

    if (eqn.psi().name().substr(0, 7) == "ILambda")
    {
        // add source term on ILambda ray
        addSupILambda(eqn);
    }
    else if (eqn.psi().name() == "Isolar")
    {
        // add source term on Isolar ray
        addSupIsolar(eqn);
    }
    else
    {
        // add source term on enthalpy equation
        addSupHe(eqn);
    }

    return;

}


void Foam::fv::absorptionEmissionSource::addSupT(fvMatrix<scalar>& eqn)
{
    DebugInformation
        << type() << ": applying source to " << eqn.psi().name() << endl;

    scalarField Ru(cells_.size(), 0);
    scalarField Rp(cells_.size(), 0);

    calculateRu(Ru);
    calculateRp(Rp);

    if (!obr_.foundObject<transportModel>("transportProperties"))
    {
        FatalErrorInFunction
            << "Could not find transportProperties object to access Cp and rho."
            << exit(FatalError);
    }

    const transportModel& transport =
        obr_.lookupObject<transportModel>("transportProperties");

    const volScalarField& Cp = transport.Cp();
    tmp<volScalarField> rho(transport.rho());

    const scalarField& V = mesh_.V();

    scalarField& TeqnSource = eqn.source();
    scalarField& TeqnDiag = eqn.diag();

    const volScalarField& T = obr_.lookupObject<volScalarField>(TName_);
    scalar TRef = basicThermo::TRefIfFound(obr_).value();

    forAll(cells_, i)
    {
        label celli = cells_[i];

        // OK!
        scalar diag =
            Rp[i]*pow3(T[celli]+TRef)/(rho()[celli]*Cp()[celli])*V[celli];
        TeqnDiag[celli] -= diag;
        TeqnSource[celli] += diag*TRef;
        // OK!
        TeqnSource[celli] -=
            Ru[i]/(rho()[celli]*Cp()[celli])*V[celli];
    }
}


void Foam::fv::absorptionEmissionSource::addSupHe
(
    fvMatrix<scalar>& eqn
)
{
    DebugInformation
        << type() << ": applying source to "
        << eqn.psi().name() << endl;

    scalarField Ru(cells_.size(), 0);
    scalarField Rp(cells_.size(), 0);

    calculateRu(Ru);
    calculateRp(Rp);

    if (!obr_.foundObject<fluidThermo>(basicThermo::dictName))
    {
        FatalErrorInFunction
            << "Could not find fluidThermo object."
            << exit(FatalError);
    }

    const fluidThermo& thermo =
        obr_.lookupObject<fluidThermo>(basicThermo::dictName);

    const scalarField& V = mesh_.V();

    scalarField& heqnSource = eqn.source();
    scalarField& heqnDiag = eqn.diag();

    tmp<volScalarField> tTAbs = thermo.TAbs();
    const volScalarField& TAbs = tTAbs();

    if (eqn.psi().name() == thermo.he().name())
    {
        const volScalarField& Cpv = thermo.Cpv();
        const volScalarField& he = eqn.psi();

        forAll(cells_, i)
        {
            label celli = cells_[i];

            heqnDiag[celli] -=
                4.0*Rp[i]*pow3(TAbs[celli])/(Cpv[celli])*V[celli];

            heqnSource[celli] -=
                (
                    Ru[i]-Rp[i]*pow3(TAbs[celli])
                   *(TAbs[celli]-4.0*he[celli]/Cpv[celli])
                )*V[celli];
        }
    }
    else
    {
        const volScalarField& T = eqn.psi();

        forAll(cells_, i)
        {
            label celli = cells_[i];

            heqnDiag[celli] -=
                4.0*Rp[i]*pow3(TAbs[celli])*V[celli];

            heqnSource[celli] -=
                (
                    Ru[i] - Rp[i]*pow3(TAbs[celli])*(TAbs[celli]-4.0*T[celli])
                )*V[celli];
        }
    }
}

void Foam::fv::absorptionEmissionSource::addSupILambda
(
    fvMatrix<scalar>& eqn
)
{
    DebugInformation
        << type() << ": applying source to " << eqn.psi().name() << endl;

    const radiationModel& radiation =
        obr_.lookupObject<radiationModel>("radiationProperties");

    const radiationModels::fvDOM& dom =
        refCast<const radiationModels::fvDOM>(radiation);

    if (!dom.participating())
    {
        WarningInFunction
            << " absorptionEmissionSource is only supported by"
            << " the participating media solver. "
            << " This will lead to unphysical results."
            << " Consider swiching to partecipating on "
            << " in radiationProperties. "<< endl;
    }

    //get rayId
    label rayId = -1;
    label lambdaId = -1;
    dom.setRayIdLambdaId(eqn.psi().name(), rayId, lambdaId);
    scalar omega_ = dom.IRay(rayId).omega();

    const scalarField& V = mesh_.V();
    const scalarField& T = obr_.lookupObject<volScalarField>(TName_);

    scalarField& IeqnSource = eqn.source();
    scalarField& IeqnDiag = eqn.diag();

    forAll(cells_, i)
    {
        const dimensionedScalar sigma_ = constant::physicoChemical::sigma;
        label celli = cells_[i];

        // absorption contribution which reduces I along ds - OK!
        IeqnDiag[celli] -= a_*omega_*V[celli];

        // emission contribution which increases I along ds - OK!
        IeqnSource[celli] -=
            1.0/constant::mathematical::pi*omega_*
            (
                  e_*sigma_.value()*pow4(T[celli]+radiation.TRef().value())
                + E_->value(obr_.time().timeOutputValue())/4
            )*V[celli];
    }

}


void Foam::fv::absorptionEmissionSource::addSupIsolar
(
    fvMatrix<scalar>& eqn
)
{
    DebugInformation
        << type() << ": applying source to " << eqn.psi().name() << endl;

    const scalarField& V = mesh_.V();

    scalarField& IseqnDiag = eqn.diag();

    forAll(cells_, i)
    {
        label celli = cells_[i];
        IseqnDiag[celli] -= a_*V[celli];
    }
}

// ************************************************************************* //
