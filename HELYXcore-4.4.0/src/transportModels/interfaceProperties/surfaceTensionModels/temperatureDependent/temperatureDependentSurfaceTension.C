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
    (c) 2017-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "surfaceTensionModels/temperatureDependent/temperatureDependentSurfaceTension.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionModels
{
    defineTypeNameAndDebug(temperatureDependent, 0);
    addToRunTimeSelectionTable
    (
        surfaceTensionModel,
        temperatureDependent,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::temperatureDependent::temperatureDependent
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    surfaceTensionModel(mesh),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    sigma_(Function1<scalar>::New("sigma", dict))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::temperatureDependent::~temperatureDependent()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::surfaceTensionModels::temperatureDependent::sigma() const
{
    tmp<volScalarField> tsigma(volScalarField::New("sigma", mesh_, dimSigma));
    volScalarField& sigma = tsigma.ref();

    basicThermo& thermo = basicThermo::lookupOrCreate(mesh_);
    tmp<volScalarField> TAbs = thermo.TAbs();

    sigma.field() = sigma_->value(TAbs().field());

    volScalarField::Boundary& sigmaBf = sigma.boundaryFieldRef();
    const volScalarField::Boundary& TBf = TAbs().boundaryField();

    forAll(sigmaBf, patchi)
    {
        sigmaBf[patchi] = sigma_->value(TBf[patchi]);
    }

    return tsigma;
}


bool Foam::surfaceTensionModels::temperatureDependent::readDict
(
    const dictionary& dict
)
{
    const dictionary& sigmaDict = surfaceTensionModel::sigmaDict(dict);

    TName_ = sigmaDict.lookupOrDefault<word>("T", "T");
    sigma_ = Function1<scalar>::New("sigma", sigmaDict);

    return true;
}


bool Foam::surfaceTensionModels::temperatureDependent::writeData
(
    Ostream& os
) const
{
    if (surfaceTensionModel::writeData(os))
    {
        os  << sigma_() << token::END_STATEMENT << nl;
        return os.good();
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
