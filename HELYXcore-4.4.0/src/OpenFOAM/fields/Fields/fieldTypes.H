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
    (c) 2019-2022 Engys Ltd.

Description
    Include the header files for all the primitive types that Fields are
    instantiated for.

\*---------------------------------------------------------------------------*/

#ifndef fieldTypes_H
#define fieldTypes_H

#include "primitives/ints/label/label.H"
#include "primitives/Scalar/scalar/scalar.H"
#include "primitives/Vector/vector/vector.H"
#include "primitives/SphericalTensor/sphericalTensor/sphericalTensor.H"
#include "primitives/SymmTensor/symmTensor/symmTensor.H"
#include "primitives/Tensor/tensor/tensor.H"
#include "primitives/Tensor4thRank/tensor4thRank/tensor4thRank.H"
#include "primitives/triad/triad.H"
#include "include/macros.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define CAPITALIZE_label Label
#define CAPITALIZE_scalar Scalar
#define CAPITALIZE_vector Vector
#define CAPITALIZE_sphericalTensor SphericalTensor
#define CAPITALIZE_symmTensor SymmTensor
#define CAPITALIZE_tensor Tensor

#define FOR_ALL_FIELD_TYPES(Macro, ...)                                        \
    Macro(scalar, __VA_ARGS__)                                                 \
    Macro(vector, __VA_ARGS__)                                                 \
    Macro(sphericalTensor, __VA_ARGS__)                                        \
    Macro(symmTensor, __VA_ARGS__)                                             \
    Macro(tensor, __VA_ARGS__)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
