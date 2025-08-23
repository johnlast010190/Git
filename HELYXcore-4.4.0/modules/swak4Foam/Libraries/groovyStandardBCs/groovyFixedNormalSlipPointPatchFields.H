/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) ICE Stroemungsfoschungs GmbH

Contributors/Copyright:
    2009-2010, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef groovyFixedNormalSlipPointPatchFields_H
#define groovyFixedNormalSlipPointPatchFields_H

#include "groovyFixedNormalSlipPointPatchField.H"
#include "fields/Fields/fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef FOAM_HAS_NO_POINT_PATCH_FIELD_TYPEDEF
    // Rebuilds the call
    // makePointPatchFieldTypedefs(groovyFixedNormalSlip);

    typedef groovyFixedNormalSlipPointPatchField<scalar> groovyFixedNormalSlipPointPatchScalarField;
    typedef groovyFixedNormalSlipPointPatchField<vector> groovyFixedNormalSlipPointPatchVectorField;
    typedef groovyFixedNormalSlipPointPatchField<sphericalTensor> groovyFixedNormalSlipPointPatchSphericalTensorField;
    typedef groovyFixedNormalSlipPointPatchField<symmTensor> groovyFixedNormalSlipPointPatchSymmTensorField;
    typedef groovyFixedNormalSlipPointPatchField<tensor> groovyFixedNormalSlipPointPatchTensorField;

#ifdef FOAM_DEV_ADDITIONAL_TENSOR_TYPES
    typedef groovyFixedNormalSlipPointPatchField<symmTensor4thOrder> groovyFixedNormalSlipPointPatchSymmTensor4thOrderField;
    typedef groovyFixedNormalSlipPointPatchField<diagTensor> groovyFixedNormalSlipPointPatchDiagTensorField;
#endif

#else
    makePointPatchFieldTypedefs(groovyFixedNormalSlip);
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
