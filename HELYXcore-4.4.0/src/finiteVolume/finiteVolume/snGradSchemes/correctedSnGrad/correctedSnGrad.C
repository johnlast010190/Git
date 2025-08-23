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
    (c) 2023 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "finiteVolume/snGradSchemes/correctedSnGrad/correctedSnGrad.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "finiteVolume/gradSchemes/gaussGrad/gaussGrad.H"

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::correctedSnGrad<Type>::~correctedSnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::fv::correctedSnGrad<Type>::fullGradCorrection
(
    const VolField<Type>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    const word& gradSchemeName
    (
        (forceGradSchemeName_ == word::null)
      ? word("grad(" + vf.name() + ')')
      : forceGradSchemeName_
    );

    // construct SurfaceField<Type>
    tmp<SurfaceField<Type>> tssf =
        linear<typename outerProduct<vector, Type>::type>(mesh).dotInterpolate
        (
            mesh.nonOrthCorrectionVectors(),
            gradScheme<Type>::New
            (
                mesh,
                vf.db(),
                mesh.schemes().gradScheme(gradSchemeName)
            )().grad(vf, gradSchemeName)
        );
    tssf.ref().rename("snGradCorr(" + vf.name() + ')');

    return tssf;
}


template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::fv::correctedSnGrad<Type>::correction
(
    const VolField<Type>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    // construct SurfaceField<Type>
    tmp<SurfaceField<Type>> tssf
    (
        SurfaceField<Type>::New
        (
            "snGradCorr(" + vf.name() + ')',
            vf.db(),
            mesh,
            vf.dimensions()*mesh.nonOrthDeltaCoeffs().dimensions()
        )
    );
    SurfaceField<Type>& ssf = tssf.ref();
    ssf.setOriented();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        ssf.replace
        (
            cmpt,
            correctedSnGrad<typename pTraits<Type>::cmptType>(mesh)
           .fullGradCorrection(vf.component(cmpt))
        );
    }

    return tssf;
}


// ************************************************************************* //
