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
    (c) 2023 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fvc/fvcSnGrad.H"
#include "fvMesh/fvMesh.H"
#include "finiteVolume/snGradSchemes/snGradScheme/snGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<SurfaceField<Type>>
snGrad
(
    const VolField<Type>& vf,
    const word& name,
    const word& forceGradSchemeName
)
{
    return fv::snGradScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().snGradScheme(name),
        forceGradSchemeName
    )().snGrad(vf);
}


template<class Type>
tmp<SurfaceField<Type>>
snGrad
(
    const tmp<VolField<Type>>& tvf,
    const word& name,
    const word& forceGradSchemeName
)
{
    tmp<SurfaceField<Type>> SnGrad
    (
        fvc::snGrad(tvf(), name, forceGradSchemeName)
    );
    tvf.clear();
    return SnGrad;
}


template<class Type>
tmp<SurfaceField<Type>>
snGrad
(
    const VolField<Type>& vf
)
{
    return fvc::snGrad(vf, "snGrad(" + vf.name() + ')');
}


template<class Type>
tmp<SurfaceField<Type>>
snGrad
(
    const tmp<VolField<Type>>& tvf
)
{
    tmp<SurfaceField<Type>> SnGrad
    (
        fvc::snGrad(tvf())
    );
    tvf.clear();
    return SnGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
