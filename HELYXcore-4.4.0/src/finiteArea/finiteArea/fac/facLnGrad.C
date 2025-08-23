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

\*---------------------------------------------------------------------------*/

#include "finiteArea/fac/facLnGrad.H"
#include "faMesh/faMesh.H"
#include "finiteArea/lnGradSchemes/lnGradScheme/lnGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fac
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
lnGrad
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name,
    const word& forceGradSchemeName
)
{
    return fa::lnGradScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().lnGradScheme(name),
        forceGradSchemeName
    )().lnGrad(vf);
}


template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
lnGrad
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
    const word& name,
    const word& forceGradSchemeName
)
{
    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> LnGrad
    (
        fac::lnGrad(tvf(), name, forceGradSchemeName)
    );
    tvf.clear();
    return LnGrad;
}


template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
lnGrad
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fac::lnGrad(vf, "lnGrad(" + vf.name() + ')');
}


template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
lnGrad
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
)
{
    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> LnGrad
    (
        fac::lnGrad(tvf())
    );
    tvf.clear();
    return LnGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fac

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
