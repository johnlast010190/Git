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

Description
    Abstract base class for surface interpolation schemes.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fv/fv.H"
#include "interpolation/surfaceInterpolation/multivariateSchemes/multivariateSurfaceInterpolationScheme/multivariateSurfaceInterpolationScheme.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::multivariateSurfaceInterpolationScheme<Type>::
multivariateSurfaceInterpolationScheme
(
    const fvMesh& mesh,
    const multivariateSurfaceInterpolationScheme<Type>::fieldTable& vtfs,
    const surfaceScalarField&,
    Istream&
)
:
    mesh_(mesh),
    fields_(vtfs)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::multivariateSurfaceInterpolationScheme<Type>>
Foam::multivariateSurfaceInterpolationScheme<Type>::New
(
    const fvMesh& mesh,
    const multivariateSurfaceInterpolationScheme<Type>::fieldTable& vtfs,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
)
{
    if (fv::debug)
    {
        InfoInFunction
            << "Constructing surfaceInterpolationScheme<Type>" << endl;
    }

    const word schemeName(schemeData);

    typename IstreamConstructorTable::iterator constructorIter =
        IstreamConstructorTable_().find(schemeName);

    if (constructorIter == IstreamConstructorTable_().end())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Unknown discretisation scheme " << schemeName << nl << nl
            << exit(FatalIOError);
    }

    return constructorIter->second(mesh, vtfs, faceFlux, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::multivariateSurfaceInterpolationScheme<Type>::
~multivariateSurfaceInterpolationScheme()
{}


// ************************************************************************* //
