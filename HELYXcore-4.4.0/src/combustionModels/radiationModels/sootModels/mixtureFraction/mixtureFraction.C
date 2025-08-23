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
    (c) 2013-2020 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "radiationModels/sootModels/mixtureFraction/mixtureFraction.H"
#include "singleStepCombustion/singleStepCombustion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace sootModels
{
    defineTypeNameAndDebug(mixtureFraction, 0);
    addToRunTimeSelectionTable(sootModel, mixtureFraction, dictionary);
}
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::sootModels::mixtureFraction::mixtureFraction
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& modelType
)
:
    sootModel(dict, mesh, modelType),
    soot_
    (
        IOobject
        (
            "soot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    coeffsDict_(dict.subOrEmptyDict(modelType + "Coeffs")),
    nuSoot_(coeffsDict_.lookup<scalar>("nuSoot")),
    Wsoot_(coeffsDict_.lookup<scalar>("Wsoot")),
    sootMax_(-1),
    mappingFieldName_
    (
        coeffsDict_.lookupOrDefault<word>("mappingField", "none")
    ),
    mapFieldMax_(1)
{
    const combustionModels::singleStepCombustion& combustion =
        mesh.lookupObject<combustionModels::singleStepCombustion>
        (
            combustionModel::combustionPropertiesName
        );

    const basicSpecieMixture& mixture = combustion.mixture();
    const reaction& singleReaction = combustion.singleReaction();

    scalar totalMol = 0;
    forAll(singleReaction.rhs(), i)
    {
        const scalar stoichCoeff = singleReaction.rhs()[i].stoichCoeff;
        totalMol += mag(stoichCoeff);
    }

    totalMol += nuSoot_;

    scalarList Xi(singleReaction.rhs().size());

    scalar Wm = 0;
    forAll(singleReaction.rhs(), i)
    {
        const label speciei = singleReaction.rhs()[i].index;
        const scalar stoichCoeff = singleReaction.rhs()[i].stoichCoeff;
        Xi[i] = mag(stoichCoeff)/totalMol;
        Wm += Xi[i]*mixture.W(speciei);
    }

    scalarList Yprod0(mixture.species().size(), 0.0);

    forAll(singleReaction.rhs(), i)
    {
        const label speciei = singleReaction.rhs()[i].index;
        Yprod0[speciei] = mixture.W(speciei)/Wm*Xi[i];
    }

    const scalar XSoot = nuSoot_/totalMol;
    Wm += XSoot*Wsoot_;

    sootMax_ = XSoot*Wsoot_/Wm;

    Info<< "Maximum soot mass concentrations: " << sootMax_ << nl;

    if (mappingFieldName_ == "none")
    {
        const label index = singleReaction.rhs()[0].index;
        mappingFieldName_ = mixture.Y(index).name();
    }

    const label mapFieldIndex = mixture.species()[mappingFieldName_];

    mapFieldMax_ = Yprod0[mapFieldIndex];
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::sootModels::mixtureFraction::~mixtureFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiationModels::sootModels::mixtureFraction::correct()
{
    const volScalarField& mapField =
        mesh_.lookupObject<volScalarField>(mappingFieldName_);

    soot_ = sootMax_*(mapField/mapFieldMax_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
