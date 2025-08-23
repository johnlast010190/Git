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
    (C) 2019 OpenCFD Ltd.
\*---------------------------------------------------------------------------*/

#ifndef polySurfacePointFieldsFwd_H
#define polySurfacePointFieldsFwd_H

#include "fields/Fields/fieldTypes.H"
#include "polySurface/polySurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
class DimensionedField;

class polySurfacePointGeoMesh;

typedef Foam::DimensionedField<label, polySurfacePointGeoMesh>
    polySurfacePointLabelField;

typedef Foam::DimensionedField<scalar, polySurfacePointGeoMesh>
    polySurfacePointScalarField;

typedef Foam::DimensionedField<vector, polySurfacePointGeoMesh>
    polySurfacePointVectorField;

typedef Foam::DimensionedField<sphericalTensor, polySurfacePointGeoMesh>
    polySurfacePointSphericalTensorField;

typedef Foam::DimensionedField<symmTensor, polySurfacePointGeoMesh>
    polySurfacePointSymmTensorField;

typedef Foam::DimensionedField<tensor, polySurfacePointGeoMesh>
    polySurfacePointTensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
