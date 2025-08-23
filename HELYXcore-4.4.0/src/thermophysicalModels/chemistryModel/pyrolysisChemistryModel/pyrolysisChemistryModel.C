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
    (c) 2013-2019 OpenFOAM Foundation
    (c) 2022-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "pyrolysisChemistryModel/pyrolysisChemistryModel.H"
#include "reaction/Reactions/Reaction/Reaction.H"
#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SolidThermo, class GasThermo>
Foam::pyrolysisChemistryModel<SolidThermo, GasThermo>::pyrolysisChemistryModel
(
    const solidMulticomponentThermo& thermo
)
:
    solidChemistryModel<SolidThermo>(thermo),
    pyrolisisGases_
    (
        dynamic_cast<const Reaction<SolidThermo>&>
        (
            this->reactions_[0]
        ).gasSpecies()
    ),
    nGases_(pyrolisisGases_.size()),
    nSpecie_(this->Ys_.size() + nGases_),
    RRg_(nGases_),
    Ys0_(this->nSolids_)
{
    // create the fields for the chemistry sources
    forAll(this->RRs_, fieldi)
    {
        IOobject header
        (
            this->Ys_[fieldi].name() + "0",
            thermo.db().time().timeName(),
            thermo.db(),
            IOobject::NO_READ
        );
        // check if field exists and can be read
        if (header.typeHeaderOk<volScalarField>(true))
        {
            Ys0_.set
            (
                fieldi,
                new volScalarField
                (
                    IOobject
                    (
                        this->Ys_[fieldi].name() + "0",
                        thermo.db().time().timeName(),
                        thermo.db(),
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    this->mesh()
                )
            );
        }
        else
        {
            volScalarField Y0Default
            (
                IOobject
                (
                    "Y0Default",
                    thermo.db().time().timeName(),
                    thermo.db(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh()
            );
            Ys0_.set
            (
                fieldi,
                new volScalarField
                (
                    IOobject
                    (
                        this->Ys_[fieldi].name() + "0",
                        thermo.db().time().timeName(),
                        thermo.db(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    Y0Default
                )
            );
            // Calculate inital values of Ysi0 = rho*delta*Yi
            Ys0_[fieldi].primitiveFieldRef() =
                this->solidThermo().rho()
               *max(this->Ys_[fieldi], scalar(0.001))*this->mesh().V();
        }
    }

    forAll(RRg_, fieldi)
    {
        RRg_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RRg." + pyrolisisGases_[fieldi],
                    thermo.db().time().timeName(),
                    thermo.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, 0)
            )
        );
    }

    Info<< "pyrolysisChemistryModel: " << nl
        << indent << "Number of solids = " << this->nSolids_ << nl
        << indent << "Number of gases = " << nGases_ << nl;
    forAll(this->reactions_, i)
    {
        Info<< dynamic_cast<const Reaction<SolidThermo>&>(this->reactions_[i])
            << nl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class SolidThermo, class GasThermo>
Foam::pyrolysisChemistryModel<SolidThermo, GasThermo>::
~pyrolysisChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SolidThermo, class GasThermo>
Foam::scalarField Foam::pyrolysisChemistryModel<SolidThermo, GasThermo>::omega
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    const bool updateC0
) const
{
    scalar omegaf, omegar;

    scalarField om(nEqns(), 0.0);

    forAll(this->reactions_, i)
    {
        const Reaction<SolidThermo>& R =
            dynamic_cast<const Reaction<SolidThermo>&>(this->reactions_[i]);

        scalar omegai = omega(R, p, T, c, li, omegaf, omegar);

        scalar rhoL = 0.0;
        scalar sr = 0.0;
        if (basicThermo::dictName == basicThermo::matDictName)
        {
            forAll(R.lhs(), s)
            {
                label si = R.lhs()[s].index;
                om[si] -= omegai;
                rhoL = this->rhoModels_[si].value(p, T);
            }
            forAll(R.rhs(), s)
            {
                label si = R.rhs()[s].index;
                scalar rhoR = this->rhoModels_[si].value(p, T);
                sr = rhoR/rhoL;
                om[si] += sr*omegai;

                if (updateC0)
                {
                    Ys0_[si][li] += sr*omegai;
                }
            }
        }
        else
        {
            forAll(R.lhs(), s)
            {
                label si = R.lhs()[s].index;
                om[si] -= omegai;
                rhoL = this->solidThermo_[si].rho(p, T);
            }
            forAll(R.rhs(), s)
            {
                label si = R.rhs()[s].index;
                scalar rhoR = this->solidThermo_[si].rho(p, T);
                sr = rhoR/rhoL;
                om[si] += sr*omegai;

                if (updateC0)
                {
                    Ys0_[si][li] += sr*omegai;
                }
            }
        }
        forAll(R.grhs(), g)
        {
            om[R.grhs()[g].index + this->nSolids_] += (1.0 - sr)*omegai;
        }
    }

    return om;
}


template<class SolidThermo, class GasThermo>
Foam::scalar Foam::pyrolysisChemistryModel<SolidThermo, GasThermo>::omega
(
    const Reaction<SolidThermo>& R,
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalar& omegaf,
    scalar& omegar
) const
{
    scalarField c1(nSpecie_, 0.0);

    for (label i=0; i<nSpecie_; i++)
    {
        c1[i] = max(0.0, c[i]);
    }

    scalar kf = R.kf(p, T, c1, li);

    const label Nl = R.lhs().size();

    for (label s=0; s<Nl; s++)
    {
        label si = R.lhs()[s].index;
        const scalar exp = R.lhs()[si].exponent;

        kf *= pow(c1[si]/Ys0_[si][li], exp)*(Ys0_[si][li]);
    }

    return kf;
}


template<class SolidThermo, class GasThermo>
void Foam::pyrolysisChemistryModel<SolidThermo, GasThermo>::derivatives
(
    const scalar time,
    const scalarField& c,
    const label li,
    scalarField& dcdt
) const
{
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    dcdt = 0.0;

    dcdt = omega(p, T, c, li);

    //Total mass concentration
    scalar cTot = 0.0;
    for (label i=0; i<this->nSolids_; i++)
    {
        cTot += c[i];
    }

    scalar newCp = 0.0;
    scalar newhi = 0.0;

    if (basicThermo::dictName == basicThermo::matDictName)
    {
        for (label i=0; i<this->nSolids_; i++)
        {
            scalar dYidt = dcdt[i]/cTot;
            scalar Yi = c[i]/cTot;
            newCp += Yi*this->CpModels_[i].value(p, T);
            newhi -= dYidt*this->hfModels_[i].value(p, T);
        }
    }
    else
    {
        for (label i=0; i<this->nSolids_; i++)
        {
            scalar dYidt = dcdt[i]/cTot;
            scalar Yi = c[i]/cTot;
            newCp += Yi*this->solidThermo_[i].Cp(p, T);
            newhi -= dYidt*this->solidThermo_[i].hf();
        }
    }

    scalar dTdt = newhi/newCp;
    scalar dtMag = min(500.0, mag(dTdt));
    dcdt[nSpecie_] = dTdt*dtMag/(mag(dTdt) + 1.0e-10);

    // dp/dt = ...
    dcdt[nSpecie_ + 1] = 0.0;
}


template<class SolidThermo, class GasThermo>
void Foam::pyrolysisChemistryModel<SolidThermo, GasThermo>::jacobian
(
    const scalar t,
    const scalarField& c,
    const label li,
    scalarField& dcdt,
    scalarSquareMatrix& dfdc
) const
{
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    scalarField c2(nSpecie_, 0.0);

    for (label i=0; i<this->nSolids_; i++)
    {
        c2[i] = max(c[i], 0.0);
    }

    for (label i=0; i<nEqns(); i++)
    {
        for (label j=0; j<nEqns(); j++)
        {
            dfdc(i, j) = 0.0;
        }
    }

    // length of the first argument must be nSolids
    dcdt = omega(p, T, c2, li);

    for (label ri=0; ri<this->reactions_.size(); ri++)
    {
        const Reaction<SolidThermo>& R = this->reactions_[ri];

        scalar kf0 = R.kf(p, T, c2, li);

        forAll(R.lhs(), j)
        {
            label sj = R.lhs()[j].index;
            scalar kf = kf0;
            forAll(R.lhs(), i)
            {
                label si = R.lhs()[i].index;
                scalar exp = R.lhs()[i].exponent;
                if (i == j)
                {
                    if (exp < 1.0)
                    {
                        kf *=
                            c2[si] > SMALL
                          ? exp*pow(c2[si] + VSMALL, exp - 1.0) : 0;
                    }
                    else
                    {
                        kf *= exp*pow(c2[si], exp - 1.0);
                    }
                }
                else
                {
                    Info<< "Solid reactions have only elements on slhs"
                        << endl;
                    kf = 0.0;
                }
            }

            forAll(R.lhs(), i)
            {
                dfdc[R.lhs()[i].index][sj] -= kf;
            }
            forAll(R.rhs(), i)
            {
                dfdc[R.rhs()[i].index][sj] += kf;
            }
        }
    }

    // calculate the dcdT elements numerically
    scalar delta = 1.0e-8;
    scalarField dcdT0 = omega(p , T - delta, c2, li);
    scalarField dcdT1 = omega(p, T + delta, c2, li);

    for (label i=0; i<nEqns(); i++)
    {
        dfdc[i][nSpecie_] = 0.5*(dcdT1[i] - dcdT0[i])/delta;
    }

}


template<class SolidThermo, class GasThermo>
Foam::label Foam::pyrolysisChemistryModel<SolidThermo, GasThermo>::nEqns() const
{
    // nEqns = number of solids + gases + temperature + pressure
    return (nSpecie_ + 2);
}


template<class SolidThermo, class GasThermo>
void Foam::pyrolysisChemistryModel<SolidThermo, GasThermo>::calculate()
{
    if (!this->chemistry_)
    {
        return;
    }

    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->solidThermo().rho()
    );

    forAll(this->RRs_, i)
    {
        this->RRs_[i].field() = 0.0;
    }

    forAll(RRg_, i)
    {
        RRg_[i].field() = 0.0;
    }

    forAll(rho, celli)
    {
        const scalar delta = this->mesh().V()[celli];

        if (this->reactingCells_[celli])
        {
            scalar rhoi = rho[celli];
            scalar Ti = this->solidThermo().T()[celli];
            scalar pi = this->solidThermo().p()[celli];

            scalarField c(nSpecie_, 0.0);
            for (label i=0; i<this->nSolids_; i++)
            {
                c[i] = rhoi*this->Ys_[i][celli]*delta;
            }

            const scalarField dcdt = omega(pi, Ti, c, celli, true);

            forAll(this->RRs_, i)
            {
                this->RRs_[i][celli] = dcdt[i]/delta;
            }

            forAll(RRg_, i)
            {
                RRg_[i][celli] = dcdt[this->nSolids_ + i]/delta;
            }
        }
    }
}


template<class SolidThermo, class GasThermo>
Foam::scalar Foam::pyrolysisChemistryModel<SolidThermo, GasThermo>::solve
(
    const scalar deltaT
)
{
    scalar deltaTMin = GREAT;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->solidThermo().rho()
    );

    forAll(this->RRs_, i)
    {
        this->RRs_[i].field() = 0.0;
    }
    forAll(RRg_, i)
    {
        RRg_[i].field() = 0.0;
    }

    const scalarField& T = this->solidThermo().T();
    const scalarField& p = this->solidThermo().p();

    scalarField c(nSpecie_, 0.0);
    scalarField c0(nSpecie_, 0.0);
    scalarField dc(nSpecie_, 0.0);
    scalarField delta(this->mesh().V());

    forAll(rho, celli)
    {
        if (this->reactingCells_[celli])
        {
            scalar rhoi = rho[celli];
            scalar pi = p[celli];
            scalar Ti = T[celli];

            for (label i=0; i<this->nSolids_; i++)
            {
                c[i] = rhoi*this->Ys_[i][celli]*delta[celli];
            }

            c0 = c;

            // Initialise time progress
            scalar timeLeft = deltaT;

            // calculate the chemical source terms
            while (timeLeft > SMALL)
            {
                scalar dt = timeLeft;
                this->solve(pi, Ti, c, celli, dt, this->deltaTChem_[celli]);
                timeLeft -= dt;
            }

            deltaTMin = min(this->deltaTChem_[celli], deltaTMin);
            dc = c - c0;

            forAll(this->RRs_, i)
            {
                this->RRs_[i][celli] = dc[i]/(deltaT*delta[celli]);
            }

            forAll(RRg_, i)
            {
                RRg_[i][celli] = dc[this->nSolids_ + i]/(deltaT*delta[celli]);
            }

            // Update Ys0_
            dc = omega(pi, Ti, c0, celli, true);
        }
    }

    // Don't allow the time-step to change more than a factor of 2
    deltaTMin = min(deltaTMin, 2*deltaT);

    return deltaTMin;
}


// ************************************************************************* //
