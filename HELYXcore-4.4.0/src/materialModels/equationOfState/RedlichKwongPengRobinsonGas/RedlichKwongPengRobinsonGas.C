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
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "RedlichKwongPengRobinsonGas.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RedlichKwongPengRobinsonGas, 0);
    addToRunTimeSelectionTable(materialModel, RedlichKwongPengRobinsonGas, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RedlichKwongPengRobinsonGas::RedlichKwongPengRobinsonGas
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
    RedlichKwongPengRobinsonGas::read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RedlichKwongPengRobinsonGas::incompressible() const
{
    return p_->isConst();
}


bool Foam::RedlichKwongPengRobinsonGas::isochoric() const
{
    return T_->isConst() && p_->isConst();
}


void Foam::RedlichKwongPengRobinsonGas::updateTable(const word& modelName)
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


Foam::baseModels<Foam::scalar>* Foam::RedlichKwongPengRobinsonGas::castScalarModel
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


bool Foam::RedlichKwongPengRobinsonGas::read()
{
    initialise("p");
    initialise("T");

    Tc_ = dict_->lookup<scalar>("Tc");
    Vc_ = dict_->lookup<scalar>("Vc");
    Pc_ = dict_->lookup<scalar>("Pc");
    omega_ = dict_->lookup<scalar>("omega");
    Zc_ = Pc_*Vc_/(RR*Tc_);

    // For calculation of EoS coefficients (see Kim et al., 2012)
    const scalar Zc = Zc_*1e-3;
    const scalar cz = 1.168;
    const scalar d1 = 0.428363;
    const scalar d2 = 18.496215;
    const scalar d3 = 0.338426;
    const scalar d4 = 0.660000;
    const scalar d5 = 789.723105;
    const scalar d6 = 2.512392;

    d1_ = (d1 + d2*exp(d4*log(d3 - cz*Zc)) + d5*exp(d6*log(d3 - cz*Zc)));
    d2_ = (1 - d1_)/(1 + d1_);

    d1PlusD2_ = d1_ + d2_;
    d1TimesD2_ = d1_*d2_;
    d1MinusD2_ = d1_ - d2_;

    const scalar dRkpr = (1 + sqr(d1_))/(1 + d1_);
    const scalar yRkpr = 1 + cbrt(2*(1 + d1_)) + cbrt(4/(1 + d1_));

    const scalar aCoeff =
        (3*sqr(yRkpr) + 3*yRkpr*dRkpr + sqr(dRkpr) + dRkpr - 1)
       /sqr(3*yRkpr + dRkpr - 1);

    const scalar bCoeff = 1/(3*yRkpr + dRkpr - 1);

    a_ = aCoeff*sqr(RR*Tc_)/Pc_;
    b_ = bCoeff*RR*Tc_/Pc_;

    const scalar a1 =  66.125;
    const scalar a0 = -23.359;
    const scalar b1 = -40.594;
    const scalar b0 =  16.855;
    const scalar c1 =  5.27345;
    const scalar c0 = -0.25826;

    const scalar coeff1 = cz*Zc*a1 + a0;
    const scalar coeff2 = cz*Zc*b1 + b0;
    const scalar coeff3 = cz*Zc*c1 + c0;

    // kappa (or c_alpha)
    // NOTE that the formumlation differs from PR and SRK
    kappa_ = coeff1*sqr(omega_) + coeff2*omega_ + coeff3;

    return true;
}


// ************************************************************************* //
