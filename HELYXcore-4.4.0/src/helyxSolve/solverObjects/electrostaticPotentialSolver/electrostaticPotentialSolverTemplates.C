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
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2019-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/constraint/empty/emptyFvPatchField.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchField.H"

template<class Type>
Foam::tmp<Foam::VolField<Type>>
Foam::fv::electrostaticPotentialSolver::createSigma
(
    autoPtr<Function1<Type>>& sigmaVsTPtr
)
{
    // This field has to be registered in the database
    // Hence New constructor is not used here
    return tmp<VolField<Type>>
    (
        new VolField<Type>
        (
            IOobject
            (
                "electrical_sigma",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            sqr(dimCurrent)/dimPower/dimLength
        )
    );
}


template<class Type>
void Foam::fv::electrostaticPotentialSolver::updateSigma
(
    const autoPtr<Function1<Type>>& sigmaVsTPtr,
    VolField<Type>& sigma
) const
{
    typedef VolField<Type> VolFieldType;

    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    // Internal field
    forAll(sigma, i)
    {
        sigma[i] = sigmaVsTPtr->value(T[i]);
    }

    // Boundary field
    typename VolFieldType::Boundary& bf = sigma.boundaryFieldRef();
    forAll(bf, patchi)
    {
        fvPatchField<Type>& pf = bf[patchi];
        if (!isA<emptyFvPatch>(pf.patch()))
        {
            const scalarField& Tbf = T.boundaryField()[patchi];
            forAll(pf, facei)
            {
                pf[facei] = sigmaVsTPtr->value(Tbf[facei]);
            }
        }
    }
}


// ************************************************************************* //
