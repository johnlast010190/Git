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
    (c) 2014-2017 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "PengRobinsonGasLaw.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PengRobinsonGasLaw, 0);
    addToRunTimeSelectionTable(materialModel, PengRobinsonGasLaw, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PengRobinsonGasLaw::PengRobinsonGasLaw
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, specieName, name)
{
    sMod_.setSize(modelsEnumSize_);
    PengRobinsonGasLaw::read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::PengRobinsonGasLaw::incompressible() const
{
    return p_->isConst();
}


bool Foam::PengRobinsonGasLaw::isochoric() const
{
    return T_->isConst() && p_->isConst();
}


void Foam::PengRobinsonGasLaw::updateTable(const word& modelName)
{
    // Model names
    const wordList modList
    ({
        rhoModel::typeName,
        hContributionModel::typeName,
        CpContributionModel::typeName,
        eContributionModel::typeName,
        CvContributionModel::typeName,
        spContributionModel::typeName,
        psiModel::typeName,
        CpMCvModel::typeName
    });

    dep_.setSize(modList.size());

    // Required dependencies
    const word& RName = RModel::typeName;
    const word& WName = WModel::typeName;
    const List<wordList> depNames
    ({
        {RName}, // rho
        {RName}, // hContribution
        {WName}, // CpContribution
        {RName}, // eContribution
        {WName}, // CvContribution
        {RName}, // spContribution
        {RName}, // psi
        {RName}  // CpMCv
    });

    // Dependency indices
    const List<labelList> depInds
    ({
        {Rind}, // rho
        {Rind}, // hContribution
        {Wind}, // CpContribution
        {Rind}, // eContribution
        {Wind}, // CvContribution
        {Rind}, // spContribution
        {Rind}, // psi
        {Rind}  // CpMCv
    });

    // Create links to dependent models
    fill(modelName, modList, depNames, depInds);
}


Foam::baseModels<Foam::scalar>* Foam::PengRobinsonGasLaw::castScalarModel
(
    const word& modelName
)
{
    castMaterial(modelName, rho)
    castMaterial(modelName, hContribution)
    castMaterial(modelName, CpContribution)
    castMaterial(modelName, eContribution)
    castMaterial(modelName, CvContribution)
    castMaterial(modelName, spContribution)
    castMaterial(modelName, svContribution)
    castMaterial(modelName, psi)
    castMaterial(modelName, Z)
    castMaterial(modelName, CpMCv)
    castMaterial(modelName, alphav)
    return nullptr;
}


bool Foam::PengRobinsonGasLaw::read()
{
    initialise("p");
    initialise("T");

    Tc_= dict_->lookup<scalar>("Tc");
    Vc_ = dict_->lookup<scalar>("Vc");
    Pc_ = dict_->lookup<scalar>("Pc");
    omega_ = dict_->lookup<scalar>("omega");
    Zc_ = Pc_*Vc_/(RR*Tc_);
    a_ = 0.45724*sqr(RR*Tc_)/Pc_;
    b_ = 0.07780*RR*Tc_/Pc_;
    kappa_ = 0.37464 + 1.54226*omega_ - 0.26992*sqr(omega_);

    return true;
}


// ************************************************************************* //
