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
    (c) 2015-2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "rhoThermo/legacyRhoThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicRhoThermo, class MixtureType>
void Foam::legacyRhoThermo<BasicRhoThermo, MixtureType>::calculate
(
    const volScalarField& p,
    volScalarField& T,
    volScalarField& he,
    volScalarField& psi,
    volScalarField& rho,
    volScalarField& mu,
    volScalarField& kappa,
    volScalarField& Cp,
    volScalarField& Cv,
    const bool doOldTimes
)
{
    // Note: update oldTimes before current time so that if T.oldTime() is
    // created from T, it starts from the unconverted T
    if (doOldTimes && (p.nOldTimes() || T.nOldTimes()))
    {
        calculate
        (
            p.oldTime(),
            T.oldTime(),
            he.oldTime(),
            psi.oldTime(),
            rho.oldTime(),
            mu.oldTime(),
            kappa.oldTime(),
            Cp.oldTime(),
            Cv.oldTime(),
            true
        );
    }

    const scalarField& hCells = he.primitiveField();
    const scalarField& pCells = p.primitiveField();

    scalarField& TCells = T.primitiveFieldRef();
    scalarField& CpCells = Cp.primitiveFieldRef();
    scalarField& CvCells = Cv.primitiveFieldRef();
    scalarField& psiCells = psi.primitiveFieldRef();
    scalarField& rhoCells = rho.primitiveFieldRef();
    scalarField& muCells = mu.primitiveFieldRef();
    scalarField& kappaCells = kappa.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoMixtureType& thermoMixture =
            this->cellThermoMixture(celli);

        const typename MixtureType::transportMixtureType& transportMixture =
            this->cellTransportMixture(celli, thermoMixture);

        const scalar p = pCells[celli] + this->pRefValue();

        TCells[celli] = thermoMixture.The(hCells[celli], p, TCells[celli]);
        CpCells[celli] = thermoMixture.Cp(p, TCells[celli]);
        CvCells[celli] = thermoMixture.Cv(p, TCells[celli]);
        psiCells[celli] = thermoMixture.psi(p, TCells[celli]);
        rhoCells[celli] = thermoMixture.rho(p, TCells[celli]);
        muCells[celli] = transportMixture.mu(p, TCells[celli]);
        kappaCells[celli] = transportMixture.kappa(p, TCells[celli]);
    }

    const volScalarField::Boundary& pBf = p.boundaryField();
    volScalarField::Boundary& TBf = T.boundaryFieldRef();
    volScalarField::Boundary& CpBf = Cp.boundaryFieldRef();
    volScalarField::Boundary& CvBf = Cv.boundaryFieldRef();
    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();
    volScalarField::Boundary& rhoBf = rho.boundaryFieldRef();
    volScalarField::Boundary& heBf = he.boundaryFieldRef();
    volScalarField::Boundary& muBf = mu.boundaryFieldRef();
    volScalarField::Boundary& kappaBf = kappa.boundaryFieldRef();

    forAll(pBf, patchi)
    {
        const fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];
        const bool fixesValue = pT.fixesValue();

        forAll(pT, facei)
        {
            const typename MixtureType::thermoMixtureType& thermoMixture =
                this->patchFaceThermoMixture(patchi, facei);
            const typename MixtureType::transportMixtureType&
                transportMixture =
                this->patchFaceTransportMixture(patchi, facei, thermoMixture);
            const scalar p = pp[facei] + this->pRefValue();
            if (fixesValue)
            {
                phe[facei] = thermoMixture.he(p, pT[facei]);
            }
            else
            {
                pT[facei] = thermoMixture.The(phe[facei], p, pT[facei]);
            }
            pCp[facei] = thermoMixture.Cp(p, pT[facei]);
            pCv[facei] = thermoMixture.Cv(p, pT[facei]);
            ppsi[facei] = thermoMixture.psi(p, pT[facei]);
            prho[facei] = thermoMixture.rho(p, pT[facei]);
            pmu[facei] = transportMixture.mu(p, pT[facei]);
            pkappa[facei] = transportMixture.kappa(p, pT[facei]);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicRhoThermo, class MixtureType>
Foam::legacyRhoThermo<BasicRhoThermo, MixtureType>::legacyRhoThermo
(
    const objectRegistry& obr,
    const word& phaseName
)
:
    legacyBasicThermo<BasicRhoThermo, MixtureType>(obr, phaseName)
{
    calculate
    (
        this->p_,
        this->T_,
        this->he_,
        this->psi_,
        this->rho_,
        this->mu_,
        this->kappa_,
        this->Cp_,
        this->Cv_,
        true                    // Create old time fields
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicRhoThermo, class MixtureType>
Foam::legacyRhoThermo<BasicRhoThermo, MixtureType>::~legacyRhoThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicRhoThermo, class MixtureType>
void Foam::legacyRhoThermo<BasicRhoThermo, MixtureType>::correct()
{
    DebugInFunction << endl;

    calculate
    (
        this->p_,
        this->T_,
        this->he_,
        this->psi_,
        this->rho_,
        this->mu_,
        this->kappa_,
        this->Cp_,
        this->Cv_,
        false           // No need to update old times
    );

    DebugInFunction << "Finished" << endl;
}


// ************************************************************************* //
