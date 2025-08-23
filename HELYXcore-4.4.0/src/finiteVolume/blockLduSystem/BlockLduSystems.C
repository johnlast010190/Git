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
    (c) 2022-2024 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "blockLduSystem/BlockLduSystem.H"
#include "blockLduSystem/BlockLduSystem.C"
#include "blockLduSystem/BlockLduSystemOperations.C"

// * * * * * * * * * * * * * * * Specialisations * * * * * * * * * * * * * * //

namespace Foam
{
template<>
void BlockLduSystem<vector, scalar>::addContinuityCoupledBC
(
    const VolField<vector>& U,
    const volScalarField& p,
    const surfaceTensorField& rDf
)
{
    const typename VolField<vector>::Boundary& bFields = U.boundaryField();
    forAll(bFields, patchi)
    {
        bFields[patchi].addContinuityCoupledBC(*this, p, rDf);
    }
}
}


// * * * * * * * * * * * * * * * Instantiations  * * * * * * * * * * * * * * //

namespace Foam
{
    template class BlockLduSystem<scalar, scalar>;
    template class BlockLduSystem<scalar, vector>;
    template class BlockLduSystem<vector, scalar>;
    template class BlockLduSystem<vector, vector>;
    template class BlockLduSystem<vector, tensor>;
    template class BlockLduSystem<tensor, tensor>;
    template class BlockLduSystem<symmTensor, symmTensor>;
    template class BlockLduSystem<sphericalTensor, sphericalTensor>;
    template class BlockLduSystem<vector1, vector1>;
    template class BlockLduSystem<vector2, vector2>;
    template class BlockLduSystem<vector4, vector4>;
    template class BlockLduSystem<vector5, vector5>;
}

// ************************************************************************* //
