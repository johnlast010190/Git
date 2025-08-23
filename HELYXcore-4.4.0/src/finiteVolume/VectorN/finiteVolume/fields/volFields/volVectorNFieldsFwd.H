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

Class
    volVectorNFields

Description

SourceFiles
    volVectorNFields.C

\*---------------------------------------------------------------------------*/

#ifndef volVectorNFieldsFwd_H
#define volVectorNFieldsFwd_H

#include "fields/volFields/volFieldsFwd.H"

#include "VectorN/Fields/VectorNFieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class volMesh;

template<class Type>
class fvPatchField;

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricField;

#define doMakeTypedefs(type, Type, args...)         \
    typedef VolField<type> vol##Type##Field;

forAllVectorNTypes(doMakeTypedefs)

forAllTensorNTypes(doMakeTypedefs)

forAllDiagTensorNTypes(doMakeTypedefs)

forAllSphericalTensorNTypes(doMakeTypedefs)

#undef doMakeTypedefs

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
