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

Class
    BlockLduMatrix

Description
    Template specialisation for symmTensor block matrix

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    symmTensorBlockLduMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef symmTensorBlockLduMatrix_H
#define symmTensorBlockLduMatrix_H

#include "matrices/blockLduMatrix/BlockLduMatrix/BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void BlockLduMatrix<symmTensor>::sumDiag();

template<>
void BlockLduMatrix<symmTensor>::negSumDiag();

template<>
void BlockLduMatrix<symmTensor>::check() const;

template<>
void BlockLduMatrix<symmTensor>::relax
(
    const symmTensorField& x,
    symmTensorField& b,
    const scalar alpha
);

template<>
void BlockLduMatrix<symmTensor>::operator*=(const scalarField& sf);

template<>
void BlockLduMatrix<symmTensor>::AmulCore
(
    symmTensorField& mul,
    const symmTensorField& x
) const;

template<>
void BlockLduMatrix<symmTensor>::TmulCore
(
    symmTensorField& mul,
    const symmTensorField& x
) const;

template<>
void BlockLduMatrix<symmTensor>::segregateB
(
    symmTensorField& mul,
    const symmTensorField& x
) const;


template<>
tmp<symmTensorField>
BlockLduMatrix<symmTensor>::H(const symmTensorField& x) const;

template<>
tmp<symmTensorField>
BlockLduMatrix<symmTensor>::faceH(const symmTensorField& x) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
