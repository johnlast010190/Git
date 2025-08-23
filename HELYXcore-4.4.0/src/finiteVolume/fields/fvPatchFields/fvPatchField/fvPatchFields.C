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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2020-2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "cfdTools/general/include/fvCFD.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchField.C"
#include "VectorN/Fields/VectorNFieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeFvPatchField(fvPatchLabelField)
makeFvPatchField(fvPatchScalarField)
makeFvPatchField(fvPatchVectorField)
makeFvPatchField(fvPatchSphericalTensorField)
makeFvPatchField(fvPatchSymmTensorField)
makeFvPatchField(fvPatchTensorField)

#define instantiateFvPatchField(Type, ...) template class fvPatchField<Type>;

FOR_ALL_FIELD_TYPES(instantiateFvPatchField)
forAllVectorNTypes(instantiateFvPatchField)
forAllTensorNTypes(instantiateFvPatchField)
forAllDiagTensorNTypes(instantiateFvPatchField)
forAllSphericalTensorNTypes(instantiateFvPatchField)


// * * * * * * * * * * * * * Static Member Functionsp * * * * * * * * * * * * //

template<class Type>
template<class Type2>
Foam::tmp<Foam::Field<Type2>>
Foam::fvPatchField<Type>::boundarySources
(
    const fvPatchField<Type2>& pf,
    const Field<Type2>& f,
    Field<Type2>& df
)
{
    fv::options& fvOptions =
        pf.db().template lookupObjectRef<fv::options>
        (
            fv::options::typeName
        );
    return
        fvOptions.boundarySources<Type>
        (
            pf.internalField().name(),
            pf.patch().index(),
            f,
            df
        );
}

// Note: Partial specialisations are implemented in fvPatchField.C


// Explicit instantiations
template Foam::tmp<Foam::Field<scalar>>
Foam::fvPatchField<scalar>::boundarySources
(
    const fvPatchField<scalar>&,
    const Field<scalar>&,
    Field<scalar>&
);

template Foam::tmp<Foam::Field<vector>>
Foam::fvPatchField<vector>::boundarySources
(
    const fvPatchField<vector>&,
    const Field<vector>&,
    Field<vector>&
);

template Foam::tmp<Foam::Field<sphericalTensor>>
Foam::fvPatchField<sphericalTensor>::boundarySources
(
    const fvPatchField<sphericalTensor>&,
    const Field<sphericalTensor>&,
    Field<sphericalTensor>&
);

template Foam::tmp<Foam::Field<symmTensor>>
Foam::fvPatchField<symmTensor>::boundarySources
(
    const fvPatchField<symmTensor>&,
    const Field<symmTensor>&,
    Field<symmTensor>&
);

template Foam::tmp<Foam::Field<tensor>>
Foam::fvPatchField<tensor>::boundarySources
(
    const fvPatchField<tensor>&,
    const Field<tensor>&,
    Field<tensor>&
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
