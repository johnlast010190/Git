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

#include "finiteVolume/fv/fv.H"
#include "finiteVolume/snGradSchemes/snGradScheme/snGradScheme.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "containers/HashTables/HashTable/HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
tmp<snGradScheme<Type>> snGradScheme<Type>::New
(
    const fvMesh& mesh,
    Istream& schemeData,
    const word& forceGradSchemeName
)
{
    if (fv::debug)
    {
        InfoInFunction << "Constructing snGradScheme<Type>" << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Discretisation scheme not specified"
            << endl << endl
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    typename MeshConstructorTable::iterator constructorIter =
        MeshConstructorTable_().find(schemeName);

    if (constructorIter == MeshConstructorTable_().end())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Unknown discretisation scheme "
            << schemeName << nl << nl
            << exit(FatalIOError);
    }

    return constructorIter->second(mesh, schemeData, forceGradSchemeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
snGradScheme<Type>::~snGradScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<SurfaceField<Type>> snGradScheme<Type>::snGrad
(
    const VolField<Type>& vf,
    const tmp<surfaceScalarField>& tdeltaCoeffs,
    const word& snGradName
)
{
    const fvMesh& mesh = vf.mesh();

    // construct SurfaceField<Type>
    tmp<SurfaceField<Type>> tsf
    (
        SurfaceField<Type>::New
        (
            snGradName + "(" + vf.name() + ')',
            vf.db(),
            vf.mesh(),
            vf.dimensions()*tdeltaCoeffs().dimensions()
        )
    );
    SurfaceField<Type>& ssf = tsf.ref();
    ssf.setOriented();

    // set reference to difference factors array
    const scalarField& deltaCoeffs = tdeltaCoeffs();

    // owner/neighbour addressing
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        ssf[facei] =
            deltaCoeffs[facei]*(vf[neighbour[facei]] - vf[owner[facei]]);
    }

    typename SurfaceField<Type>::Boundary& ssfbf = ssf.boundaryFieldRef();

    forAll(vf.boundaryField(), patchi)
    {
        const fvPatchField<Type>& pvf = vf.boundaryField()[patchi];

        if (pvf.coupled())
        {
            ssfbf[patchi] = pvf.snGrad(tdeltaCoeffs().boundaryField()[patchi]);
        }
        else
        {
            ssfbf[patchi] = pvf.snGrad();
        }
    }

    return tsf;
}


template<class Type>
tmp<SurfaceField<Type>> snGradScheme<Type>::sndGrad
(
    const VolField<Type>& vf,
    const word& sndGradName
)
{
    return snGrad(vf, vf.mesh().nonOrthDeltaCoeffs(), sndGradName);
}


template<class Type>
tmp<SurfaceField<Type>> snGradScheme<Type>::snGrad
(
    const VolField<Type>& vf
) const
{
    tmp<SurfaceField<Type>> tsf(snGrad(vf, deltaCoeffs(vf)));

    if (corrected())
    {
        tsf.ref() += correction(vf);
    }
    return tsf;
}


template<class Type>
tmp<SurfaceField<Type>> snGradScheme<Type>::snGrad
(
    const tmp<VolField<Type>>& tvf
) const
{
    tmp<SurfaceField<Type>> tsf(snGrad(tvf()));
    tvf.clear();
    return tsf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
