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
    (c) 2016-2017 Wikki Ltd

\*---------------------------------------------------------------------------*/

#include "faMesh/faMesh.H"
#include "fields/areaFields/areaFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTemplate2TypeNameAndDebug(areaScalarField::Internal, 0);
defineTemplate2TypeNameAndDebug(areaVectorField::Internal, 0);
defineTemplate2TypeNameAndDebug(areaSphericalTensorField::Internal, 0);
defineTemplate2TypeNameAndDebug(areaSymmTensorField::Internal, 0);
defineTemplate2TypeNameAndDebug(areaTensorField::Internal, 0);

defineTemplateTypeNameAndDebug(areaScalarField, 0);
defineTemplateTypeNameAndDebug(areaVectorField, 0);
defineTemplateTypeNameAndDebug(areaSphericalTensorField, 0);
defineTemplateTypeNameAndDebug(areaSymmTensorField, 0);
defineTemplateTypeNameAndDebug(areaTensorField, 0);

template<>
tmp<GeometricField<scalar, faPatchField, areaMesh>>
GeometricField<scalar, faPatchField, areaMesh>::component
(
    const direction
) const
{
    return *this;
}

template<>
void GeometricField<scalar, faPatchField, areaMesh>::replace
(
    const direction,
    const GeometricField<scalar, faPatchField, areaMesh>& gsf
)
{
    this->forceAssign(gsf);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
