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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "Gulder.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{
    defineTypeNameAndDebug(Gulder, 0);
    addToRunTimeSelectionTable(XiEqModel, Gulder, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiEqModels::Gulder::Gulder
(
    const dictionary& XiEqProperties,
    const psiuMulticomponentThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su
)
:
    XiEqModel(XiEqProperties, thermo, turbulence, Su),
    XiEqCoef_(XiEqModelCoeffs_.lookup<scalar>("XiEqCoef")),
    SuMin_(0.01*Su.average()),
    uPrimeCoef_(XiEqModelCoeffs_.lookup<scalar>("uPrimeCoef")),
    subGridSchelkin_(XiEqModelCoeffs_.lookup<bool>("subGridSchelkin"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiEqModels::Gulder::~Gulder()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiEqModels::Gulder::XiEq() const
{
    volScalarField up(sqrt((2.0/3.0)*turbulence_.k()));
    tmp<volScalarField> tepsilon(turbulence_.epsilon());
    const volScalarField& epsilon = tepsilon();

    if (subGridSchelkin_)
    {
        up.primitiveFieldRef() += calculateSchelkinEffect(uPrimeCoef_);
    }

    volScalarField tauEta(sqrt(mag(thermo_.muu()/(thermo_.rhou()*epsilon))));

    volScalarField Reta
    (
        up
      / (
            sqrt(epsilon*tauEta)
          + dimensionedScalar(up.dimensions(), 1e-8)
        )
    );

    return (1.0 + XiEqCoef_*sqrt(up/(Su_ + SuMin_))*Reta);
}


bool Foam::XiEqModels::Gulder::read(const dictionary& XiEqProperties)
{
    XiEqModel::read(XiEqProperties);

    XiEqCoef_ = XiEqModelCoeffs_.lookup<scalar>("XiEqCoef");
    uPrimeCoef_ = XiEqModelCoeffs_.lookup<scalar>("uPrimeCoef");
    subGridSchelkin_ = XiEqModelCoeffs_.lookup<bool>("subGridSchelkin");

    return true;
}


// ************************************************************************* //
