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
    Typedefs for block matrices

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef blockLduMatrices_H
#define blockLduMatrices_H

#include "matrices/blockLduMatrix/BlockLduMatrix/scalarBlockLduMatrix.H"
#include "matrices/blockLduMatrix/BlockLduMatrix/sphericalTensorBlockLduMatrix.H"
#include "matrices/blockLduMatrix/BlockLduMatrix/symmTensorBlockLduMatrix.H"
#include "matrices/blockLduMatrix/BlockLduMatrix/tensorBlockLduMatrix.H"
#include "containers/HashTables/Map/Map.H"
#include "matrices/blockLduMatrix/BlockLduMatrix/BlockConstraint/blockConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef BlockLduMatrix<scalar> blockScalarMatrix;
typedef BlockLduMatrix<vector> blockVectorMatrix;
typedef BlockLduMatrix<sphericalTensor> blockSphericalTensorMatrix;
typedef BlockLduMatrix<symmTensor> blockSymmTensorMatrix;
typedef BlockLduMatrix<tensor> blockTensorMatrix;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
