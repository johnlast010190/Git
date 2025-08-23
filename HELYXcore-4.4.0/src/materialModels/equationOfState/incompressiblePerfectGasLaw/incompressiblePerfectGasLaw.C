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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "incompressiblePerfectGasLaw.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressiblePerfectGasLaw, 0);
    addToRunTimeSelectionTable(materialModel, incompressiblePerfectGasLaw, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressiblePerfectGasLaw::incompressiblePerfectGasLaw
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
    incompressiblePerfectGasLaw::read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::incompressiblePerfectGasLaw::incompressible() const
{
    return true;
}


bool Foam::incompressiblePerfectGasLaw::isochoric() const
{
    return T_->isConst();
}


void Foam::incompressiblePerfectGasLaw::updateTable(const word& modelName)
{
    // Model names
    const wordList modList
    ({
        rhoModel::typeName,
        hContributionModel::typeName,
        CpContributionModel::typeName,
        ZModel::typeName
    });
    dep_.setSize(modList.size());

    // Required dependencies
    const wordList depNames(modList.size(), RModel::typeName);

    // Dependency indices
    const labelList depInds(modList.size(), RInd);

    // Create up links to dependent models
    fill(modelName, modList, depNames, depInds);
}


Foam::baseModels<Foam::scalar>*
Foam::incompressiblePerfectGasLaw::castScalarModel
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


bool Foam::incompressiblePerfectGasLaw::read()
{
    initialise("T");
    initialise("p");

    // For const p => pRef = PConst + Pref
    if (dict_->found("pRef"))
    {
        if (dict_->found("pConst"))
        {
            FatalErrorInFunction
                << "Please specify 'pRef' (absolute pressure) "
                << "or 'pConst' (gauge pressure) but not both"
                << exit(FatalError);
        }
        else
        {
            DeprecationIOWarningInFunction
            (
                dict_,
                "pRef",
                "entry",
                40200,
                "Please specify gauge pressure with 'pConst' instead"
            );
            pRef_ = dict_->lookup<scalar>("pRef");
        }
    }
    else if (p_->isConst() && !dict_->found("pConst"))
    {
        pRef_ = p_->constant().value() + p_->offset().value();
    }
    else
    {
        pRef_ = dict_->lookup<scalar>("pConst") + p_->offset().value();
    }
    Info<< nl << nl << "Reference pressure in incompressible perfect gas law: "
        << pRef_ << " Pa" << nl;
    return true;
}


// ************************************************************************* //
