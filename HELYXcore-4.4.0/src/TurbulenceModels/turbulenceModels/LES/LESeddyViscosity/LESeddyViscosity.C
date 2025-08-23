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
    (c) 2016 OpenCFD Ltd.
    (c) 2013-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "LESeddyViscosity.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
LESeddyViscosity<BasicTurbulenceModel>::LESeddyViscosity
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    eddyViscosity<LESModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    Ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ck",
            this->coeffDict_,
            0.094
        )
    ),
    Ce_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ce",
            this->coeffDict_,
            1.048
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool LESeddyViscosity<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<LESModel<BasicTurbulenceModel>>::read())
    {
        Ck_.readIfPresent(this->coeffDict());
        Ce_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField> LESeddyViscosity<BasicTurbulenceModel>::epsilon() const
{
    tmp<volScalarField> tk(this->k());


    tmp<volScalarField> tepsilon
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("epsilon", this->U_.group()),
                this->runTime_.timeName(),
                this->db_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            Ce_*tk()*sqrt(tk())/this->delta(),
            zeroGradientFvPatchField<scalar>::typeName
        )
    );
    volScalarField& epsilon = tepsilon.ref();
    epsilon.correctBoundaryConditions();

    return tepsilon;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> LESeddyViscosity<BasicTurbulenceModel>::omega() const
{
    tmp<volScalarField> tk(this->k());
    tmp<volScalarField> tomega
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("omega", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->db_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            (Ce_/Ck_)*sqrt(tk())/this->delta(),
            zeroGradientFvPatchField<scalar>::typeName
        )
    );
    volScalarField& omega = tomega.ref();
    omega.correctBoundaryConditions();

    return tomega;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
