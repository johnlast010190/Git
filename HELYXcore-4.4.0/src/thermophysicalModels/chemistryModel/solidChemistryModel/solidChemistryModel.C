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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "solidChemistryModel/solidChemistryModel.H"
#include "mixtures/multicomponentMixture/multicomponentMixture.H"
#include "multiphaseThermo/multiphaseThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SolidThermo>
Foam::solidChemistryModel<SolidThermo>::solidChemistryModel
(
    const solidMulticomponentThermo& thermo
)
:
    basicSolidChemistryModel(thermo),
    ODESystem(),
    Ys_
    (
        const_cast<PtrList<volScalarField>&>
        (
            this->solidThermo().composition().Y()
        )
    ),
    solidThermo_
    (
        dynamic_cast<const multicomponentMixture<SolidThermo>&>
        (
            this->solidThermo()
        ).specieThermos()
    ),
    reactions_
    (
        dynamic_cast<const multicomponentMixture<SolidThermo>&>
        (
            this->solidThermo()
        ).species(),
        solidThermo_,
        this->mesh(),
        *this
    ),
    nSolids_(Ys_.size()),
    nReaction_(reactions_.size()),
    RRs_(nSolids_),
    reactingCells_(this->mesh().nCells(), true)
{
    if (basicThermo::dictName == basicThermo::matDictName)
    {
        hfModels_.resize(nSolids_);
        hsModels_.resize(nSolids_);
        rhoModels_.resize(nSolids_);
        CpModels_.resize(nSolids_);
        const word phaseName(Ys_.first().group());
        forAll(Ys_, speciei)
        {
            const word specieName(Ys_[speciei].member());
            matScalarTable& speciesModels =
                thermo.materials().sTable(phaseName, specieName);
            hfModels_.set(speciei, speciesModels[hfModel::typeName]);
            hsModels_.set(speciei, speciesModels[hsModel::typeName]);
            rhoModels_.set(speciei, speciesModels[rhoModel::typeName]);
            CpModels_.set(speciei, speciesModels[CpModel::typeName]);
        }
    }

    // create the fields for the chemistry sources
    forAll(RRs_, fieldi)
    {
        RRs_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RRs." + Ys_[fieldi].name(),
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class SolidThermo>
Foam::solidChemistryModel<SolidThermo>::~solidChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SolidThermo>
Foam::tmp<Foam::volScalarField>
Foam::solidChemistryModel<SolidThermo>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        volScalarField::New
        (
            "Qdot",
            this->mesh_,
            dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
        )
    );

    if (this->chemistry_)
    {
        scalarField& Qdot = tQdot.ref();

        forAll(Ys_, i)
        {
            forAll(Qdot, celli)
            {
                const scalar hfi
                (
                    (basicThermo::dictName == basicThermo::matDictName)
                  ? hfModels_[i][celli]
                  : solidThermo_[i].hf()
                );
                Qdot[celli] -= hfi*RRs_[i][celli];
            }
        }
    }

    return tQdot;
}


template<class SolidThermo>
void Foam::solidChemistryModel<SolidThermo>::setCellReacting
(
    const label celli,
    const bool active
)
{
    reactingCells_[celli] = active;
}

// ************************************************************************* //
