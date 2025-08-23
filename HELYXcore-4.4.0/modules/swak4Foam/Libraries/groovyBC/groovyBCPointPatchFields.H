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
    (c) 1991-2008 OpenCFD Ltd.

Contributors/Copyright:
    2010, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef groovyBCPointPatchFields_H
#define groovyBCPointPatchFields_H

#include "groovyBCPointPatchField.H"
#include "fields/Fields/fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef FOAM_HAS_NO_POINT_PATCH_FIELD_TYPEDEF
    // Rebuilds the call
    // makePointPatchFieldTypedefs(groovyBC);

    typedef groovyBCPointPatchField<scalar> groovyBCPointPatchScalarField;
    typedef groovyBCPointPatchField<vector> groovyBCPointPatchVectorField;
    typedef groovyBCPointPatchField<sphericalTensor> groovyBCPointPatchSphericalTensorField;
    typedef groovyBCPointPatchField<symmTensor> groovyBCPointPatchSymmTensorField;
    typedef groovyBCPointPatchField<tensor> groovyBCPointPatchTensorField;

#ifdef FOAM_DEV_ADDITIONAL_TENSOR_TYPES
    typedef groovyBCPointPatchField<symmTensor4thOrder> groovyBCPointPatchSymmTensor4thOrderField;
    typedef groovyBCPointPatchField<diagTensor> groovyBCPointPatchDiagTensorField;
#endif

#else
makePointPatchFieldTypedefs(groovyBC);
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
