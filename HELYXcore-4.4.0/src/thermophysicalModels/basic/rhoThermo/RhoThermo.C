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
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2021-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "RhoThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicRhoThermo, class MixtureType>
void Foam::RhoThermo<BasicRhoThermo, MixtureType>::calculate
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
    DebugInFunction << "Updating material models." << endl;

    // Note: update oldTimes before current time so that if T.oldTime() is
    // created from T, it starts from the unconverted T
    if (doOldTimes && (p.nOldTimes() || T.nOldTimes()))
    {
        this->materials_.updateScalarField
        (
            BasicRhoThermo::phasePropertyName("p"),
            p.oldTime()
        );
        this->materials_.updateScalarField
        (
            BasicRhoThermo::phasePropertyName("T"),
            T.oldTime()
        );
        this->materials_.updateScalarField
        (
            BasicRhoThermo::phasePropertyName(he.name()),
            he.oldTime()
        );
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

    // Update temperature internal field
    if (this->calculatesTFromhe())
    {
        if
        (
            this->materials_.isMultiphase()
         && this->phaseName_ == word::null
         && this->materials_.isPhasicVariable("T")
        )
        {
            // If T is being solved phasic, then he comes from a mixture
            // and so needs to be saved. This should probably be done in a
            // more general way through the models or maybe even the solver
            // - or global T should be turned off entirely
            he.primitiveFieldRef() =
                this->materials_
                (
                    heModel::typeName,
                    this->phaseName_
                ).primitiveField();
        }
        T.primitiveFieldRef() =
            this->materials_
            (
                TModel::typeName,
                this->phaseName_
            ).primitiveField();
    }
    else
    {
        he.primitiveFieldRef() =
            this->materials_
            (
                heModel::typeName,
                this->phaseName_
            ).primitiveField();
    }

    // Is temperature constant?
    const bool isTConst =
        obr_.found("TRef")
     && obr_.lookupObject<refScalarField>("TRef").isConst();

    // Update he and T boundary fields
    volScalarField::Boundary& TBf = T.boundaryFieldRef();
    volScalarField::Boundary& heBf = he.boundaryFieldRef();
    forAll(TBf, patchi)
    {
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];

        if (this->calculatesTFromhe())
        {
            if (pT.fixesValue() && !isTConst)
            {
                phe.forceAssign
                (
                    this->materials_
                    (
                        heModel::typeName,
                        this->phaseName_
                    ).boundaryField()[patchi]
                );
            }
            else
            {
                if
                (
                    this->materials_.isMultiphase()
                 && this->phaseName_ == word::null
                 && this->materials_.isPhasicVariable("T")
                )
                {
                    // If T is being solved phasic, then he comes from a mixture
                    // and so needs to be saved. This should probably be done in
                    // a more general way through the models or maybe even the
                    // solver - or global T should be turned off entirely
                    phe.forceAssign
                    (
                        this->materials_
                        (
                            heModel::typeName,
                            this->phaseName_
                        ).boundaryField()[patchi]
                    );
                }

                pT.forceAssign
                (
                    this->materials_
                    (
                        TModel::typeName,
                        this->phaseName_
                    ).boundaryField()[patchi]
                );
            }
        }
        else
        {
            phe.forceAssign
            (
                this->materials_
                (
                    heModel::typeName,
                    this->phaseName_
                ).boundaryField()[patchi]
            );
        }
    }
    Cp.forceAssign(this->materials_(CpModel::typeName, this->phaseName_)());
    Cv.forceAssign(this->materials_(CvModel::typeName, this->phaseName_)());
    rho.forceAssign(this->materials_(rhoModel::typeName, this->phaseName_)());
    psi.forceAssign(this->materials_(psiModel::typeName, this->phaseName_)());

    // In case of a phase mixture of solid and liquid, ignore nonexistent mu
    if
    (
        this->phaseName_ == word::null
     || this->materials_.sTable
        (
            this->phaseName_,
            word::null
        ).found(muModel::typeName)
    )
    {
        mu.forceAssign(this->materials_(muModel::typeName, this->phaseName_)());
    }
    kappa.forceAssign
    (
        this->materials_(kappaModel::typeName, this->phaseName_)()
    );

    // Update back when old times is true and we are doing current timestep
    // which isn't 0 (nOldTimes() = 0 means there were no old times
    // and update isn't necessary)
    if
    (
        doOldTimes
     && this->p_.nOldTimes()
     && (this->p_.nOldTimes() == p.nOldTimes())
    )
    {
        this->materials_.updateScalarField
        (
            BasicRhoThermo::phasePropertyName("p"),
            this->p_
        );
        this->materials_.updateScalarField
        (
            BasicRhoThermo::phasePropertyName("T"),
            this->T_
        );
        this->materials_.updateScalarField
        (
            BasicRhoThermo::phasePropertyName(he.name()),
            this->he_
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicRhoThermo, class MixtureType>
Foam::RhoThermo<BasicRhoThermo, MixtureType>::RhoThermo
(
    const objectRegistry& obr,
    const word& phaseName
)
:
    BasicThermo<BasicRhoThermo, MixtureType>(obr, phaseName),
    obr_(obr)
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
Foam::RhoThermo<BasicRhoThermo, MixtureType>::~RhoThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicRhoThermo, class MixtureType>
void Foam::RhoThermo<BasicRhoThermo, MixtureType>::correct()
{
    BasicRhoThermo::correct();

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
}


// ************************************************************************* //
