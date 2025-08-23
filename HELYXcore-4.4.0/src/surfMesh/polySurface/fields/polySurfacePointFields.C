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

#include "polySurface/fields/polySurfacePointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const word polySurfacePointLabelField::typeName
("polySurfacePointLabelField");

template<>
const word polySurfacePointScalarField::typeName
("polySurfacePointScalarField");

template<>
const word polySurfacePointVectorField::typeName
("polySurfacePointVectorField");

template<>
const word polySurfacePointSphericalTensorField::typeName
("polySurfacePointSphericalTensorField");

template<>
const word polySurfacePointSymmTensorField::typeName
("polySurfacePointSymmTensorField");

template<>
const word polySurfacePointTensorField::typeName
("polySurfacePointTensorField");


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
