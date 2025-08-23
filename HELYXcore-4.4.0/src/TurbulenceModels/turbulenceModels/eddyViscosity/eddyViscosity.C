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
    (c) 2013-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "eddyViscosity.H"
#include "finiteVolume/fvc/fvc.H"
#include "finiteVolume/fvm/fvm.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::eddyViscosity<BasicTurbulenceModel>::eddyViscosity
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
    linearViscousStress<BasicTurbulenceModel>
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
    nut_
    (
        IOobject
        (
            IOobject::groupName("nut", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->db_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool Foam::eddyViscosity<BasicTurbulenceModel>::read()
{
    return BasicTurbulenceModel::read();
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::eddyViscosity<BasicTurbulenceModel>::R() const
{
    tmp<volScalarField> tk(k());

    // Get list of patchField type names from k
    wordList patchFieldTypes(tk().boundaryField().types());

    // For k patchField types which do not have an equivalent for symmTensor
    // set to calculated
    forAll(patchFieldTypes, i)
    {
        if (fvPatchField<symmTensor>::patchConstructorTable_().count(patchFieldTypes[i]) == 0)
        {
            patchFieldTypes[i] = calculatedFvPatchField<symmTensor>::typeName;
        }
    }

    return tmp<volSymmTensorField>
    (
        volSymmTensorField::New
        (
            IOobject::groupName("R", this->alphaRhoPhi_.group()),
            ((2.0/3.0)*I)*tk() - (nut_)*dev(twoSymm(fvc::grad(this->U_))),
            patchFieldTypes
        )
    );
}


template<class BasicTurbulenceModel>
void Foam::eddyViscosity<BasicTurbulenceModel>::validate()
{
    correctNut();
}


template<class BasicTurbulenceModel>
void Foam::eddyViscosity<BasicTurbulenceModel>::correct()
{
    BasicTurbulenceModel::correct();
}


// ************************************************************************* //
