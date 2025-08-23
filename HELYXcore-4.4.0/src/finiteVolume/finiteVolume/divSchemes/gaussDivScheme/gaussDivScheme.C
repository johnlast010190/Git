
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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2019-2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/divSchemes/gaussDivScheme/gaussDivScheme.H"
#include "finiteVolume/fvc/fvcSurfaceIntegrate.H"
#include "fvMatrices/fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolField<typename innerProduct<vector, Type>::type>>
gaussDivScheme<Type>::fvcDiv(const VolField<Type>& vf)
{
    tmp<VolField<typename innerProduct<vector, Type>::type>> tDiv
    (
        fvc::surfaceIntegrate
        (
            this->tinterpScheme_().dotInterpolate(this->mesh_.Sf(), vf)
        )
    );

    tDiv.ref().rename("div(" + vf.name() + ')');

    return tDiv;
}


template<class Type>
tmp<BlockLduSystem<vector, typename innerProduct<vector, Type>::type>>
gaussDivScheme<Type>::fvmDiv(const VolField<Type>& vf) const
{
    FatalErrorInFunction
        << "Only scalar Type is supported"
        << exit(FatalError);

    typedef typename innerProduct<vector, Type>::type DivType;

    tmp<BlockLduSystem<vector, DivType>> tbs
    (
        new BlockLduSystem<vector, DivType>(vf.mesh())
    );

    return tbs;
}


template<class Type>
tmp<BlockLduSystem<vector, typename innerProduct<vector, Type>::type>>
gaussDivScheme<Type>::fvmBCDiv
(
    const SurfaceField<scalar>& rhof,
    const VolField<Type>& vf
) const
{
    FatalErrorInFunction
        << "Only scalar Type is supported"
        << exit(FatalError);

    typedef typename innerProduct<vector, Type>::type DivType;

    tmp<BlockLduSystem<vector, DivType>> tbs
    (
        new BlockLduSystem<vector, DivType>(vf.mesh())
    );

    return tbs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
