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
    (c) 2018-2025 Engys Ltd.
    (c) 2011-2022 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "porousAbsorptionEmissionSource.H"
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
    defineTypeNameAndDebug(porousAbsorptionEmissionSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        porousAbsorptionEmissionSource,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::porousAbsorptionEmissionSource::porousAbsorptionEmissionSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    absorptionEmissionSource(name, modelType, dict, obr),
    Tsolid_
    (
        IOobject
        (
            "Tsolid_"+cellSetName_,
            obr_.time().timeName(),
            obr_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        obr_.lookupObject<volScalarField>(TName_)
    ),
    rhoSolid_(coeffs_.lookup<scalar>("rhoSolid")),
    CpSolid_(coeffs_.lookup<scalar>("CpSolid")),
    hCoeff_(coeffs_.lookup<scalar>("hCoefficient"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::porousAbsorptionEmissionSource::addSupHe
(
    fvMatrix<scalar>& eqn
)
{
    DebugInformation
        << type() << ": applying source to " << eqn.psi().name() << endl;

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

    const fluidThermo& thermo
        (obr_.lookupObject<fluidThermo>(basicThermo::dictName));

    const scalarField& V = mesh_.V();
    const volScalarField& T = obr_.lookupObject<volScalarField>(TName_);
    scalar TRef = basicThermo::TRefIfFound(obr_).value();

    scalarField& heqnSource = eqn.source();
    scalarField& heqnDiag = eqn.diag();

    if (eqn.psi().name() == thermo.he().name())
    {
        const volScalarField& Cpv = thermo.Cpv();
        const volScalarField& he = eqn.psi();

        forAll(cells_, i)
        {
            label celli = cells_[i];

            heqnDiag[celli] -=
                4.0*Rp[i]*pow3(Tsolid_[celli]+TRef)/(Cpv[celli])*V[celli];

            heqnSource[celli] -=
                (
                    // Ru -Rp*T^4: absorbed - emitted heat flux
                    // 4*Rp*T^3*h/Cpv terms cancel out, added for stability
                    Ru[i]-Rp[i]*pow3(Tsolid_[celli]+TRef)
                    *(Tsolid_[celli]+TRef -4.0*he[celli]/Cpv[celli])
                    // heat conduction term; TRef cancel out
                    + hCoeff_*(Tsolid_[celli] - T[celli])
                )*V[celli];
        }

        //- update Tsolid_ of porous medium
        calculateTsolid();
    }
    else
    {
        forAll(cells_, i)
        {
            label celli = cells_[i];

            heqnDiag[celli] -=
                4.0*Rp[i]*pow3(Tsolid_[celli]+TRef)*V[celli];

            heqnSource[celli] -=
                (
                    // Ru -Rp*T^4; 4*T terms cancel out, added for stability
                    Ru[i] - Rp[i]*pow3(Tsolid_[celli]+TRef)*
                    (Tsolid_[celli]+TRef -4.0*(Tsolid_[celli]+TRef))
                    // heat conduction term; TRef cancel out
                    + hCoeff_*(Tsolid_[celli] - T[celli])
                )*V[celli];
        }

        calculateTsolid();
    }
}

void Foam::fv::porousAbsorptionEmissionSource::addSupILambda
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
            << " porousAbsorptionEmissionSource is only supported by"
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

    scalarField& IeqnSource = eqn.source();
    scalarField& IeqnDiag = eqn.diag();

    forAll(cells_, i)
    {
        const dimensionedScalar sigma_ = constant::physicoChemical::sigma;
        label celli = cells_[i];

        // absorption contribution which reduces I along ds
        IeqnDiag[celli] -= a_*omega_*V[celli];

        // emission contribution which increases I along ds
        IeqnSource[celli] -=
            1.0/constant::mathematical::pi*omega_*
            (
                e_*sigma_.value()*pow4(Tsolid_[celli]+radiation.TRef().value())
                + E_->value(obr_.time().timeOutputValue())/4
            )*V[celli];
    }

}


void Foam::fv::porousAbsorptionEmissionSource::calculateTsolid()
{
    const volScalarField& T = obr_.lookupObject<volScalarField>(TName_);
    const volScalarField& G_ = obr_.lookupObject<volScalarField>("G");

    label zoneIndex = mesh_.cellZones().findZoneID(cellSetName_);
    const labelList& zoneCells
    (
        mesh_.cellZones()[zoneIndex]
    );

    forAll(zoneCells, ci)
    {
        label celli = zoneCells[ci];
        // convected heat from fluid region
        scalar Qconv = hCoeff_*(T[celli] - Tsolid_[celli]);
        // absorbed heat via radiation from fluid - emitted from porous medium
        scalar Qrad = a_*G_[celli]
            -4*e_*constant::physicoChemical::sigma.value()*pow4(Tsolid_[celli]);
        Tsolid_[celli] += (mesh_.time().deltaT().value()/(rhoSolid_*CpSolid_))
            *(Qconv + Qrad);

        DebugInformation
            << "Qabsorbed " << a_*G_[celli] << ", Qemitted "
            << 4*e_*constant::physicoChemical::sigma.value()*pow4(Tsolid_[celli])
            << ", Qconvection " << Qconv << endl;
    }
}

// ************************************************************************* //
