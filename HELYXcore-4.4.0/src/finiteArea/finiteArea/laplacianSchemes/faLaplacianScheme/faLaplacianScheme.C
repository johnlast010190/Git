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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016-2024 Engys Ltd.

Description
    Abstract base class for finite area calculus laplacian schemes.

\*---------------------------------------------------------------------------*/

#include "finiteArea/fa/fa.H"
#include "containers/HashTables/HashTable/HashTable.H"
#include "interpolation/edgeInterpolation/schemes/linear/linearEdgeInterpolation.H"
#include "faMatrices/faMatrix/faMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<laplacianScheme<Type, GType>> laplacianScheme<Type, GType>::New
(
    const faMesh& mesh,
    Istream& schemeData,
    const word& forceGradSchemeName
)
{
    if (fa::debug)
    {
        Info<< "laplacianScheme<Type>::New(const faMesh&, Istream&) : "
               "constructing laplacianScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Laplacian scheme not specified" << endl << endl
            << exit(FatalIOError);
    }

    word schemeName(schemeData);

    const auto ctor =
        ctorTableLookup
        (
            "laplacian scheme", IstreamConstructorTable_(), schemeName
        );
    return ctor(mesh, schemeData, forceGradSchemeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, class GType>
laplacianScheme<Type, GType>::~laplacianScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<faMatrix<Type>>
laplacianScheme<Type, GType>::famLaplacian
(
//    const areaScalarField& gamma,
    const GeometricField<GType, faPatchField, areaMesh>& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return famLaplacian(tinterpGammaScheme_().interpolate(gamma)(), vf);
}


template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacianScheme<Type, GType>::facLaplacian
(
    const GeometricField<GType, faPatchField, areaMesh>& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return facLaplacian(tinterpGammaScheme_().interpolate(gamma)(), vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
