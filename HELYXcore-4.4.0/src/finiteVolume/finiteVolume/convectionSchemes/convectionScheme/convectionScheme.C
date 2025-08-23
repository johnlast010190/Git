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
    (c) 2019-2024 Engys Ltd.

Description
    Abstract base class for finite volume calculus convection schemes.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fv/fv.H"
#include "containers/HashTables/HashTable/HashTable.H"
#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
convectionScheme<Type>::convectionScheme(const convectionScheme& cs)
:
    tmp<convectionScheme<Type>>::refCount(),
    mesh_(cs.mesh_)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
tmp<convectionScheme<Type>> convectionScheme<Type>::New
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
)
{
    if (fv::debug)
    {
        InfoInFunction << "Constructing convectionScheme<Type>" << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Convection scheme not specified" << endl << endl
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    const auto ctor =
        ctorTableLookup
        (
            "convection scheme",
            IstreamConstructorTable_(),
            schemeName
        );
    return ctor(mesh, faceFlux, schemeData);
}


template<class Type>
tmp<convectionScheme<Type>> convectionScheme<Type>::New
(
    const fvMesh& mesh,
    const typename multivariateSurfaceInterpolationScheme<Type>::
        fieldTable& fields,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
)
{
    if (fv::debug)
    {
        InfoInFunction << "Constructing convectionScheme<Type>" << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Convection scheme not specified" << endl << endl
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    const auto ctor =
        ctorTableLookup
        (
            "convection scheme",
            MultivariateConstructorTable_(),
            schemeName
        );
    return ctor(mesh, fields, faceFlux, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
convectionScheme<Type>::~convectionScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<fvBlockMatrix<Type>> convectionScheme<Type>::fvmBDiv
(
    const surfaceScalarField& phi,
    const VolField<Type>& vf
) const
{
    FatalErrorInFunction
        << "Only gaussConvectionScheme is supported"
        << exit(FatalError);

    tmp<fvBlockMatrix<Type>> tfvbm
    (
        new fvBlockMatrix<Type>
        (
            const_cast<VolField<Type>&>(vf)
        )
    );

    return tfvbm;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void convectionScheme<Type>::operator=(const convectionScheme<Type>& cs)
{
    if (this == &cs)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
