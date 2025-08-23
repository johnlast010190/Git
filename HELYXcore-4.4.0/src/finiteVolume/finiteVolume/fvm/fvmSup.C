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

#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Su
(
    const DimensionedField<Type, volMesh>& su,
    const VolField<Type>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            dimVol*su.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    fvm.source() -= mesh.V()*su.field();

    return tfvm;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Su
(
    const tmp<DimensionedField<Type, volMesh>>& tsu,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm = fvm::Su(tsu(), vf);
    tsu.clear();
    return tfvm;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Su
(
    const tmp<VolField<Type>>& tsu,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm = fvm::Su(tsu(), vf);
    tsu.clear();
    return tfvm;
}


template<class Type>
Foam::zeroField
Foam::fvm::Su
(
    const zero&,
    const VolField<Type>& vf
)
{
    return zeroField();
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Sp
(
    const volScalarField::Internal& sp,
    const VolField<Type>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            dimVol*sp.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    fvm.diag() += mesh.V()*sp.field();

    return tfvm;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Sp
(
    const tmp<volScalarField::Internal>& tsp,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm = fvm::Sp(tsp(), vf);
    tsp.clear();
    return tfvm;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Sp
(
    const tmp<volScalarField>& tsp,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm = fvm::Sp(tsp(), vf);
    tsp.clear();
    return tfvm;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Sp
(
    const dimensionedScalar& sp,
    const VolField<Type>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            dimVol*sp.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    fvm.diag() += mesh.V()*sp.value();

    return tfvm;
}


template<class Type>
Foam::zeroField
Foam::fvm::Sp
(
    const zero&,
    const VolField<Type>&
)
{
    return zeroField();
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>>
Foam::fvm::bSp
(
    const tmp
    <
        GeometricField
        <
            typename outerProduct<vector, Type>::type,
            fvPatchField,
            volMesh
        >
    >& tgradf,
    const VolField<Type>& vf
)
{
    tmp<fvBlockMatrix<Type>> tbEq
    (
        new fvBlockMatrix<Type>
        (
            const_cast<VolField<Type>&>(vf),
            tgradf().dimensions()*dimVolume
        )
    );

    fvBlockMatrix<Type>& bs = tbEq.ref();

    typename CoeffField<Type>::squareTypeField& d = bs.diag().asSquare();

    d = tgradf()*vf.mesh().V();

    return tbEq;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>>
Foam::fvm::bSuSp
(
    const tmp
    <
        GeometricField
        <
            typename outerProduct<vector, Type>::type,
            fvPatchField,
            volMesh
        >
    >& tgradf,
    const VolField<Type>& vf
)
{
    tmp<fvBlockMatrix<Type>> tbEq
    (
        new fvBlockMatrix<Type>
        (
            const_cast<VolField<Type>&>(vf),
            tgradf().dimensions()*dimVolume
        )
    );

    fvBlockMatrix<Type>& bs = tbEq.ref();

    typename CoeffField<Type>::squareTypeField& d = bs.diag().asSquare();

    tensorField& bDiag = bs.diag().asSquare();
    vectorField& source = bs.source();

    d = tgradf()*vf.mesh().V();

    label nCoeffs = pTraits<vector>::nComponents;
    forAll(bDiag, cI)
    {
        for (label iC=0; iC<nCoeffs; iC++)
        {
            label tDiag = iC*nCoeffs+iC;
            label ofDiag1 = (nCoeffs-1)*iC + 1 + max(0, (iC-1));
            label ofDiag2 = nCoeffs*iC + 2 - max(0, (iC-1));
            if (tgradf()[cI][tDiag]<0)
            {
                source[cI].component(iC) -= (tgradf()[cI].component(tDiag)*vf.primitiveField()[cI][iC])*vf.mesh().V()[cI];
                source[cI].component(iC) -= (tgradf()[cI].component(ofDiag1)*vf.primitiveField()[cI][ofDiag1%3])*vf.mesh().V()[cI];
                source[cI].component(iC) -= (tgradf()[cI].component(ofDiag2)*vf.primitiveField()[cI][ofDiag2%3])*vf.mesh().V()[cI];
                d[cI].component(tDiag) = 0;
                d[cI].component(ofDiag1) = 0;
                d[cI].component(ofDiag2) = 0;
            }
        }
    }

    return tbEq;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::SuSp
(
    const volScalarField::Internal& susp,
    const VolField<Type>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            dimVol*susp.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    fvm.diag() += mesh.V()*max(susp.field(), scalar(0));

    fvm.source() -= mesh.V()*min(susp.field(), scalar(0))
        *vf.primitiveField();

    return tfvm;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::SuSp
(
    const tmp<volScalarField::Internal>& tsusp,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm = fvm::SuSp(tsusp(), vf);
    tsusp.clear();
    return tfvm;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::SuSp
(
    const tmp<volScalarField>& tsusp,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm = fvm::SuSp(tsusp(), vf);
    tsusp.clear();
    return tfvm;
}


template<class Type>
Foam::zeroField
Foam::fvm::SuSp
(
    const zero&,
    const VolField<Type>& vf
)
{
    return zeroField();
}


// ************************************************************************* //
