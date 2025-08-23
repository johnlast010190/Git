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

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fvc/fvcCurl.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolField<Type>> curl(const VolField<Type>& vf)
{
    word nameCurlVf = "curl(" + vf.name() + ')';

    // Gausses theorem curl
    // tmp<VolField<Type>> tcurlVf =
    //     fvc::surfaceIntegrate(vf.mesh().Sf() ^ fvc::interpolate(vf));

    // Calculate curl as the Hodge dual of the skew-symmetric part of grad
    tmp<VolField<Type>> tcurlVf =
        2.0*(*skew(fvc::grad(vf, nameCurlVf)));

    tcurlVf.ref().rename(nameCurlVf);

    return tcurlVf;
}


template<class Type>
tmp<VolField<Type>> curl(const tmp<VolField<Type>>& tvf)
{
    tmp<VolField<Type>> Curl(fvc::curl(tvf()));
    tvf.clear();
    return Curl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
