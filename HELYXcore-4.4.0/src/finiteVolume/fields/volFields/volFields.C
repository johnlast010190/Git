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
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTemplate2TypeNameAndDebug(volLabelField::Internal, 0);
defineTemplate2TypeNameAndDebug(volScalarField::Internal, 0);
defineTemplate2TypeNameAndDebug(volVectorField::Internal, 0);
defineTemplate2TypeNameAndDebug
(
    volSphericalTensorField::Internal,
    0
);
defineTemplate2TypeNameAndDebug
(
    volSymmTensorField::Internal,
    0
);
defineTemplate2TypeNameAndDebug(volTensorField::Internal, 0);

defineTemplateTypeNameAndDebug(volLabelField, 0);
defineTemplateTypeNameAndDebug(volScalarField, 0);
defineTemplateTypeNameAndDebug(volVectorField, 0);
defineTemplateTypeNameAndDebug(volSphericalTensorField, 0);
defineTemplateTypeNameAndDebug(volSymmTensorField, 0);
defineTemplateTypeNameAndDebug(volTensorField, 0);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specialization for scalar fields
template<>
tmp<VolField<scalar>>
VolField<scalar>::component
(
    const direction
) const
{
    return *this;
}


// Specialization for scalar fields
template<>
void VolField<scalar>::replace
(
    const direction,
    const VolField<scalar>& gsf
)
{
    this->forceAssign(gsf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
