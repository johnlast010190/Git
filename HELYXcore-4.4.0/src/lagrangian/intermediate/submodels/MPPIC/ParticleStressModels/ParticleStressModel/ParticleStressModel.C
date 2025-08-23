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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "submodels/MPPIC/ParticleStressModels/ParticleStressModel/ParticleStressModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ParticleStressModel, 0);
    defineRunTimeSelectionTable(ParticleStressModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ParticleStressModel::ParticleStressModel
(
    const dictionary& dict
)
:
    alphaPacked_(dict.lookup<scalar>("alphaPacked"))
{
}


Foam::ParticleStressModel::ParticleStressModel
(
    const ParticleStressModel& cm
)
:
    alphaPacked_(cm.alphaPacked_)
{
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::ParticleStressModel> Foam::ParticleStressModel::New
(
    const dictionary& dict
)
{
    word modelType(dict.lookup("type"));

    Info<< "Selecting particle stress model " << modelType << endl;

    const auto ctor =
        ctorTableLookup
        (
            "particle stress model type",
            dictionaryConstructorTable_(),
            modelType
        );
    return autoPtr<ParticleStressModel>(ctor(dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ParticleStressModel::~ParticleStressModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::ParticleStressModel::alphaPacked() const
{
    return alphaPacked_;
}


Foam::tmp<Foam::FieldField<Foam::Field, Foam::scalar>>
Foam::ParticleStressModel::tau
(
    const FieldField<Field, scalar>& alpha,
    const FieldField<Field, scalar>& rho,
    const FieldField<Field, scalar>& uRms
) const
{
    tmp<FieldField<Field, scalar>> value
    (
        new FieldField<Field, scalar>(alpha.size())
    );

    forAll(alpha, i)
    {
        value->set(i, tau(alpha[i], rho[i], uRms[i]));
    }

    return value;
}


// ************************************************************************* //
