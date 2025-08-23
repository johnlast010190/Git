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
    (c) 2004-6 H. Jasak All rights reserved

\*---------------------------------------------------------------------------*/

#include "matrices/blockLduMatrix/BlockLduMatrix/sphericalTensorBlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void Foam::BlockLduMatrix<Foam::sphericalTensor>::sumDiag()
{
    // Decoupled version
    this->decoupledSumDiag();
}


template<>
void Foam::BlockLduMatrix<Foam::sphericalTensor>::negSumDiag()
{
    // Decoupled version
    this->decoupledNegSumDiag();
}


template<>
void Foam::BlockLduMatrix<sphericalTensor>::check() const
{
    // Decoupled version
    this->decoupledCheck();
}


template<>
void Foam::BlockLduMatrix<Foam::sphericalTensor>::relax
(
    const sphericalTensorField& x,
    sphericalTensorField& b,
    const scalar alpha
)
{
    // Decoupled version
    this->decoupledRelax(x, b, alpha);
}


template<>
void Foam::BlockLduMatrix<Foam::sphericalTensor>::operator*=
(
    const scalarField& sf
)
{
    // Decoupled version
    this->decoupledMultEqOp(sf);
}


template<>
void Foam::BlockLduMatrix<Foam::sphericalTensor>::AmulCore
(
    sphericalTensorField& Ax,
    const sphericalTensorField& x
) const
{
    decoupledAmulCore(Ax, x);
}


template<>
void Foam::BlockLduMatrix<Foam::sphericalTensor>::TmulCore
(
    sphericalTensorField& Tx,
    const sphericalTensorField& x
) const
{
    // Decoupled version
    decoupledTmulCore(Tx, x);
}


template<>
void Foam::BlockLduMatrix<Foam::sphericalTensor>::segregateB
(
    sphericalTensorField&,
    const sphericalTensorField&
) const
{
    FatalErrorInFunction
        << "Requested decoupling of sphericalTensor matrix - never coupled"
        << abort(FatalError);
}


template<>
Foam::tmp<Foam::sphericalTensorField>
Foam::BlockLduMatrix<Foam::sphericalTensor>::H
(
    const sphericalTensorField& x
) const
{
    // Decoupled version
    return decoupledH(x);
}


template<>
Foam::tmp<Foam::sphericalTensorField>
Foam::BlockLduMatrix<Foam::sphericalTensor>::faceH
(
    const sphericalTensorField& x
) const
{
    // Decoupled version
    return decoupledFaceH(x);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
