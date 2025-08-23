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
    (c) 2004-2011 H. Jasak

Description
    Global functions for expansion and contraction of tensor field
    to diagonal type

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "VectorN/expandContract/ExpandTensorNField.H"
#include "VectorN/expandContract/ExpandTensorN.H"
#include "fields/Fields/Field/Field.H"
#include "VectorN/Fields/VectorNFieldTypes.H"
#include "fields/Fields/Field/FieldM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(typeF1, typeF2, FUNC)                                 \
                                                                             \
void FUNC(Field<typeF1>& f1, const UList<typeF2>& f2)                        \
{                                                                            \
    checkFields(f1, f2, #FUNC "(f1,f2)");                                    \
                                                                             \
    List_ACCESS(typeF1, f1, f1P);                                            \
    List_CONST_ACCESS(typeF2, f2, f2P);                                      \
                                                                             \
    List_FOR_ALL(f1,i)                                                       \
        FUNC(List_ELEM(f1, f1P, i), List_ELEM(f2, f2P, i));                  \
    List_END_FOR_ALL                                                         \
}                                                                            \
                                                                             \
void FUNC(Field<typeF1>& f1, const tmp<Field<typeF2>>& tf2)                 \
{                                                                            \
     FUNC(f1,tf2());                                                         \
     tf2.clear();                                                            \
}

#define ExpandFieldFunctions(tensorType, diagTensorType, sphericalTensorType, \
        vectorType, cmptType, args...)                                        \
                                                                              \
UNARY_FUNCTION(cmptType, tensorType, contractScalar)                          \
UNARY_FUNCTION(cmptType, diagTensorType, contractScalar)                      \
UNARY_FUNCTION(cmptType, sphericalTensorType, contractScalar)                 \
UNARY_FUNCTION(cmptType, vectorType, contractScalar)                          \
                                                                              \
UNARY_FUNCTION(vectorType, tensorType, contractLinear)                        \
UNARY_FUNCTION(vectorType, diagTensorType, contractLinear)                    \
UNARY_FUNCTION(vectorType, sphericalTensorType, contractLinear)               \
                                                                              \
UNARY_FUNCTION(vectorType, cmptType, expandScalar)                            \
UNARY_FUNCTION(tensorType, cmptType, expandScalar)                            \
UNARY_FUNCTION(diagTensorType, cmptType, expandScalar)                        \
UNARY_FUNCTION(sphericalTensorType, cmptType, expandScalar)                   \
                                                                              \
UNARY_FUNCTION(tensorType, vectorType, expandLinear)                          \
UNARY_FUNCTION(diagTensorType, vectorType, expandLinear)                      \
UNARY_FUNCTION(sphericalTensorType, vectorType, expandLinear)

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

forAllVectorTensorNTypes(ExpandFieldFunctions)

// template<class Cmpt, int length>
// void contractScalar
// (
//     Field<Cmpt>& res,
//     const UList<TensorN<Cmpt, length>>& f
// )
// {
//     forAll(res, i)
//     {
//         contractScalar(res[i], f[i]);
//     }
// }
//
//
// template <class Cmpt, int length>
// void contractLinear
// (
//     Field<VectorN<Cmpt, length>>& res,
//     const UList<TensorN<Cmpt, length>>& f
// )
// {
//     forAll(res, i)
//     {
//         contractLinear(res[i], f[i]);
//     }
// }
//
//
// template <class Cmpt, int length>
// void expandScalar
// (
//     Field<VectorN<Cmpt, length>>& res,
//     const UList<Cmpt>& f
// )
// {
//     forAll(res, i)
//     {
//         expandScalar(res[i], f[i]);
//     }
// }
//
//
// template <class Cmpt, int length>
// void expandScalar
// (
//     Field<TensorN<Cmpt, length>>& res,
//     const UList<Cmpt>& f
// )
// {
//     forAll(res, i)
//     {
//         expandScalar(res[i], f[i]);
//     }
// }
//
//
// template <class Cmpt, int length>
// void expandLinear
// (
//     Field<TensorN<Cmpt, length>>& res,
//     const UList<VectorN<Cmpt, length>>& f
// )
// {
//     forAll(res, i)
//     {
//         expandLinear(res[i], f[i]);
//     }
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef UNARY_FUNCTION
#undef ExpandFieldFunctions

// ************************************************************************* //
