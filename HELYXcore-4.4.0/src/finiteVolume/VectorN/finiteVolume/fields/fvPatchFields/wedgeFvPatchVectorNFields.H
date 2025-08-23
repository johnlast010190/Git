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
    (c) 2010 Ivor Clifford

\*---------------------------------------------------------------------------*/

#ifndef wedgeFvPatchVectorNFields_H
#define wedgeFvPatchVectorNFields_H

#include "VectorN/finiteVolume/fields/volFields/volVectorNFields.H"
#include "VectorN/finiteVolume/fields/fvPatchFields/fvPatchVectorNFields.H"

#include "VectorN/finiteVolume/fields/fvPatchFields/wedgeFvPatchVectorNFieldsFwd.H"
#include "fields/fvPatchFields/constraint/wedge/wedgeFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define doMakeVectorTensorNWedgeFuncDefs(type, Type, args...)       \
template<>                                                          \
tmp<Field<type >> wedgeFvPatchField<type >::snGrad() const;        \
                                                                    \
template<>                                                          \
void wedgeFvPatchField<type >::evaluate(                            \
    const Pstream::commsTypes commsType                             \
);                                                                  \
                                                                    \
template<>                                                          \
tmp<Field<type >> wedgeFvPatchField<type >::snGradTransformDiag()  \
const;

forAllVectorNTypes(doMakeVectorTensorNWedgeFuncDefs)

forAllTensorNTypes(doMakeVectorTensorNWedgeFuncDefs)

forAllDiagTensorNTypes(doMakeVectorTensorNWedgeFuncDefs)

forAllSphericalTensorNTypes(doMakeVectorTensorNWedgeFuncDefs)

#undef doMakeVectorTensorNWedgeFuncDefs

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
