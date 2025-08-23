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
    (c) 2011 OpenFOAM Foundation

Description
    Specialisation of Field\<T\> for diagTensor.

\*---------------------------------------------------------------------------*/

#include "fields/Fields/diagTensorField/diagTensorField.H"

#define TEMPLATE
#include "fields/Fields/Field/FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

UNARY_FUNCTION(diagTensor, tensor, diag)
UNARY_FUNCTION(scalar, diagTensor, tr)
UNARY_FUNCTION(sphericalTensor, diagTensor, sph)
UNARY_FUNCTION(scalar, diagTensor, det)
UNARY_FUNCTION(diagTensor, diagTensor, inv)


BINARY_OPERATOR(tensor, diagTensor, tensor, +, add)
BINARY_OPERATOR(tensor, diagTensor, tensor, -, subtract)

BINARY_TYPE_OPERATOR(tensor, diagTensor, tensor, +, add)
BINARY_TYPE_OPERATOR(tensor, diagTensor, tensor, -, subtract)

BINARY_OPERATOR(vector, vector, diagTensor, /, divide)
BINARY_TYPE_OPERATOR(vector, vector, diagTensor, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fields/Fields/Field/undefFieldFunctionsM.H"

// ************************************************************************* //
