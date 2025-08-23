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
    (c) 2012-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fvMatrices/fvMatrices.H"
#include "fields/volFields/volFields.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"

Foam::tmp<Foam::DimensionedField<Foam::vector, Foam::volMesh>>
Foam::coalCloudList::UTrans() const
{
    tmp<volVectorField::Internal> tfld
    (
        volVectorField::Internal::New
        (
            "UTransEff",
            mesh_,
            dimensionedVector(dimMass*dimVelocity, Zero)
        )
    );

    volVectorField::Internal& fld = tfld.ref();

    forAll(*this, i)
    {
        fld += operator[](i).UTrans();
    }

    return tfld;
}


Foam::tmp<Foam::fvVectorMatrix> Foam::coalCloudList::SU
(
    volVectorField& U
) const
{
    tmp<fvVectorMatrix> tfvm(new fvVectorMatrix(U, dimForce));
    fvVectorMatrix& fvm = tfvm.ref();

    forAll(*this, i)
    {
        fvm += operator[](i).SU(U);
    }

    return tfvm;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::coalCloudList::hsTrans() const
{
    tmp<volScalarField::Internal> tfld
    (
        volScalarField::Internal::New
        (
            "hsTransEff",
            mesh_,
            dimensionedScalar(dimEnergy, 0)
        )
    );

    volScalarField::Internal& fld = tfld.ref();

    forAll(*this, i)
    {
        fld += operator[](i).hsTrans();
    }

    return tfld;
}


Foam::tmp<Foam::fvScalarMatrix> Foam::coalCloudList::Sh
(
    volScalarField& hs
) const
{
    tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(hs, dimEnergy/dimTime));
    fvScalarMatrix& fvm = tfvm.ref();

    forAll(*this, i)
    {
        fvm += operator[](i).Sh(hs);
    }

    return tfvm;
}


Foam::tmp<Foam::fvScalarMatrix> Foam::coalCloudList::SYi
(
    const label ii,
    volScalarField& Yi
) const
{
    tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(Yi, dimMass/dimTime));
    fvScalarMatrix& fvm = tfvm.ref();

    forAll(*this, i)
    {
        fvm += operator[](i).SYi(ii, Yi);
    }

    return tfvm;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::coalCloudList::rhoTrans() const
{
    tmp<volScalarField::Internal> tfld
    (
        volScalarField::Internal::New
        (
            "rhoTransEff",
            mesh_,
            dimensionedScalar(dimMass, 0)
        )
    );

    volScalarField::Internal& fld = tfld.ref();

    forAll(*this, i)
    {
        forAll(operator[](i).rhoTrans(), j)
        {
            fld += operator[](i).rhoTrans()[j];
        }
    }

    return tfld;
}




Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::coalCloudList::Srho() const
{
    tmp<volScalarField::Internal> tfld
    (
        volScalarField::Internal::New
        (
            "rhoTransEff",
            mesh_,
            dimensionedScalar(dimDensity/dimTime, 0)
        )
    );

    volScalarField::Internal& fld = tfld.ref();

    forAll(*this, i)
    {
        fld += operator[](i).Srho();
    }

    return tfld;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::coalCloudList::Srho
(
    const label i
) const
{
    tmp<volScalarField::Internal> tfld
    (
        volScalarField::Internal::New
        (
            "rhoTransEff",
            mesh_,
            dimensionedScalar(dimDensity/dimTime, 0)
        )
    );

    volScalarField::Internal& fld = tfld.ref();

    forAll(*this, j)
    {
        fld += operator[](j).Srho(i);
    }

    return tfld;
}


Foam::tmp<Foam::fvScalarMatrix> Foam::coalCloudList::Srho
(
    volScalarField& rho
) const
{
    tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(rho, dimMass/dimTime));
    fvScalarMatrix& fvm = tfvm.ref();

    forAll(*this, i)
    {
        fvm += operator[](i).Srho(rho);
    }

    return tfvm;
}


// ************************************************************************* //
