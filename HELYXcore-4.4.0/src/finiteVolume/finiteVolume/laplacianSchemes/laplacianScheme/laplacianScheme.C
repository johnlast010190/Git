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
    (c) 2021 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fv/fv.H"
#include "containers/HashTables/HashTable/HashTable.H"
#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<laplacianScheme<Type, GType>> laplacianScheme<Type, GType>::New
(
    const fvMesh& mesh,
    const objectRegistry& db,
    Istream& schemeData,
    const word& forceGradSchemeName
)
{
    if (fv::debug)
    {
        InfoInFunction << "Constructing laplacianScheme<Type, GType>" << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Laplacian scheme not specified" << endl << endl
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    const auto ctor =
        ctorTableLookup
        (
            "laplacian scheme", IstreamConstructorTable_(), schemeName
        );
    return ctor(mesh, db, schemeData, forceGradSchemeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, class GType>
laplacianScheme<Type, GType>::~laplacianScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<fvMatrix<Type>>
laplacianScheme<Type, GType>::fvmLaplacian
(
    const VolField<GType>& gamma,
    const VolField<Type>& vf
)
{
    return fvmLaplacian(tinterpGammaScheme_().interpolate(gamma)(), vf);
}


template<class Type, class GType>
tmp<fvBlockMatrix<Type>>
laplacianScheme<Type, GType>::fvmBLaplacian
(
    const VolField<GType>& gamma,
    const VolField<Type>& vf
)
{
    return fvmBLaplacian(tinterpGammaScheme_().interpolate(gamma)(), vf);
}


template<class Type, class GType>
tmp<VolField<Type>>
laplacianScheme<Type, GType>::fvcLaplacian
(
    const VolField<GType>& gamma,
    const VolField<Type>& vf
)
{
    return fvcLaplacian(tinterpGammaScheme_().interpolate(gamma)(), vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
