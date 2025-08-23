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
    (c) 2017-2020 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "liquidPropertiesSurfaceTension.H"
#include "liquidProperties/liquidMixtureProperties/liquidMixtureProperties.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "basicThermo/basicThermo.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionModels
{
    defineTypeNameAndDebug(liquidPropertiesSurfaceTension, 0);
    addToRunTimeSelectionTable
    (
        surfaceTensionModel,
        liquidPropertiesSurfaceTension,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::liquidPropertiesSurfaceTension::
liquidPropertiesSurfaceTension
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    surfaceTensionModel(mesh),
    phaseName_(dict.lookup("phase"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::liquidPropertiesSurfaceTension::
~liquidPropertiesSurfaceTension()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::surfaceTensionModels::liquidPropertiesSurfaceTension::sigma() const
{
    IOdictionary thermoDict
    (
        IOobject
        (
            IOobject::groupName(basicThermo::dictName, phaseName_),
            db().time().constant(),
            db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    const liquidProperties thermo(thermoDict);

    tmp<volScalarField> tsigma(volScalarField::New("sigma", mesh_, dimSigma));
    volScalarField& sigma = tsigma.ref();

    const volScalarField& T = db().lookupObject<volScalarField>("T");
    const volScalarField& p = db().lookupObject<volScalarField>("p");

    volScalarField::Internal& sigmai = sigma;
    const volScalarField::Internal& pi = p;
    const volScalarField::Internal& Ti = T;

    forAll(sigmai, celli)
    {
        sigmai[celli] = thermo.sigma(pi[celli], Ti[celli]);
    }

    volScalarField::Boundary& sigmaBf = sigma.boundaryFieldRef();
    const volScalarField::Boundary& pBf = p.boundaryField();
    const volScalarField::Boundary& TBf = T.boundaryField();

    forAll(sigmaBf, patchi)
    {
        scalarField& sigmaPf = sigmaBf[patchi];
        const scalarField& pPf = pBf[patchi];
        const scalarField& TPf = TBf[patchi];

        forAll(sigmaPf, facei)
        {
            sigmaPf[facei] = thermo.sigma(pPf[facei], TPf[facei]);
        }
    }

    return tsigma;
}


bool Foam::surfaceTensionModels::liquidPropertiesSurfaceTension::readDict
(
    const dictionary& dict
)
{
    return true;
}


bool Foam::surfaceTensionModels::liquidPropertiesSurfaceTension::writeData
(
    Ostream& os
) const
{
    if (surfaceTensionModel::writeData(os))
    {
        return os.good();
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
