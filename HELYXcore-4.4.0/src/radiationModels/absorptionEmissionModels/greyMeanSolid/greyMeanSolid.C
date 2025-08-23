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
    (c) 2015 OpenCFD Ltd.
    (c) 2011-2020 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "absorptionEmissionModels/greyMeanSolid/greyMeanSolid.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/unitConversion/unitConversion.H"
#include "fields/fvPatchFields/basic/extrapolatedCalculated/extrapolatedCalculatedFvPatchFields.H"
#include "multiphaseThermo/multiphaseThermo.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{
    defineTypeNameAndDebug(greyMeanSolid, 0);

    addToRunTimeSelectionTable
    (
        absorptionEmissionModel,
        greyMeanSolid,
        dictionary
    );
    // Backward compatibility
    addNamedToRunTimeSelectionTable
    (
        absorptionEmissionModel,
        greyMeanSolid,
        dictionary,
        greyMeanSolidAbsorptionEmission
    );
}
}
}


// * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::radiationModels::absorptionEmissionModels::
greyMeanSolid::X(const word specie) const
{
    tmp<volScalarField> TAbs = thermo_.TAbs();
    const volScalarField& T = TAbs();
    const volScalarField& p = thermo_.p();

    tmp<scalarField> tXj(new scalarField(T.primitiveField().size(), 0.0));
    scalarField& Xj = tXj.ref();

    tmp<scalarField> tRhoInv(new scalarField(T.primitiveField().size(), 0.0));
    scalarField& rhoInv = tRhoInv.ref();

    forAll(mixture_.Y(), specieI)
    {
        const scalarField& Yi = mixture_.Y()[specieI];
        if (basicThermo::dictName == basicThermo::matDictName)
        {
            forAll(rhoInv, i)
            {
                rhoInv[i] +=
                    Yi[i]
                   /thermo_.materials()
                    (
                        rhoModel::typeName,
                        thermo_.phaseName(),
                        IOobject::member(mixture_.Y()[specieI].name())
                    )[i];
            }
        }
        else
        {
            forAll(rhoInv, iCell)
            {
                rhoInv[iCell] +=
                    Yi[iCell]/mixture_.rho(specieI, p[iCell], T[iCell]);
            }
        }
    }
    const scalarField& Yj = mixture_.Y(specie);
    const label mySpecieI = mixture_.species()[specie];
    if (basicThermo::dictName == basicThermo::matDictName)
    {
        forAll(Xj, iCell)
        {
            Xj[iCell] =
                Yj[iCell]
               /thermo_.materials()
                (
                    rhoModel::typeName,
                    thermo_.phaseName(),
                    IOobject::member(mixture_.Y()[mySpecieI].name())
                )[iCell];
        }
    }
    else
    {
        forAll(Xj, iCell)
        {
            Xj[iCell] = Yj[iCell]/mixture_.rho(mySpecieI, p[iCell], T[iCell]);
        }
    }

    return (Xj/rhoInv);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

const Foam::solidThermo&
Foam::radiationModels::absorptionEmissionModels::greyMeanSolid::thermo
(
    const fvMesh& mesh
)
{
    if
    (
        isA<multiphaseThermo>
        (
            mesh.lookupObject<basicThermo>(basicThermo::dictName)
        )
    )
    {
        return
            mesh.lookupObject<multiphaseThermo>
            (
                basicThermo::dictName
            ).sThermo();
    }
    return mesh.lookupObject<solidThermo>(basicThermo::dictName);
}


Foam::radiationModels::absorptionEmissionModels::greyMeanSolid::greyMeanSolid
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& defaultTypeName
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_((dict.optionalSubDict(defaultTypeName + "Coeffs"))),
    thermo_(thermo(mesh)),
    speciesNames_(0),
    mixture_(dynamic_cast<const basicSpecieMixture&>(thermo_)),
    solidData_(mixture_.Y().size())
{
    if (!isA<basicSpecieMixture>(thermo_))
    {
        FatalErrorInFunction
            << "Model requires a multi-component thermo package"
            << abort(FatalError);
    }

    label nFunc = 0;

    forAllConstIter(dictionary, coeffsDict_, iter)
    {
        // safety:
        if (!iter().isDict())
        {
            continue;
        }
        const word& key = iter().keyword();
        if (!mixture_.contains(key))
        {
            WarningInFunction
                << " specie: " << key << " is not found in the solid mixture"
                << nl
                << " specie is the mixture are:" << mixture_.species() << nl
                << nl << endl;
        }
        speciesNames_.insert(key, nFunc);
        const dictionary& dict = iter().dict();
        solidData_[nFunc][absorptivity] = dict.lookup<scalar>("absorptivity");
        solidData_[nFunc][emissivity] = dict.lookup<scalar>("emissivity");

        nFunc++;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::greyMeanSolid::
~greyMeanSolid()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::greyMeanSolid::calc
(
    const label propertyId
) const
{
    tmp<volScalarField> ta
    (
        volScalarField::New
        (
            "a",
            mesh(),
            dimensionedScalar(dimless/dimLength, 0),
            extrapolatedCalculatedFvPatchVectorField::typeName
        )
    );

    scalarField& a = ta.ref().primitiveFieldRef();

    forAllConstIter(HashTable<label>, speciesNames_, iter)
    {
        if (mixture_.contains(iter.key()))
        {
            a += solidData_[iter()][propertyId]*X(iter.key());
        }
    }

    ta.ref().correctBoundaryConditions();
    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::greyMeanSolid::eCont
(
    const label bandI
) const
{
   return calc(emissivity);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::greyMeanSolid::aCont
(
    const label bandI
) const
{
   return calc(absorptivity);
}

// ************************************************************************* //
