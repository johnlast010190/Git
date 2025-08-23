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
    2010-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceMesh/surfaceMesh.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchField.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchField.H"
#include "interpolations/primitivePatchInterpolation/primitivePatchInterpolation.H"
#include "interpolation/interpolation/interpolation/interpolation.H"

#include <cassert>

#include "CloudProxy.H"

namespace Foam {

    template <typename Type>
    tmp<Field<Type >> CloudValueExpressionDriver::getFluidField(
        word const&name
    )
    {
        Field<vector> positions(proxy_->getPositions());
        Field<label> cells(proxy_->getCells());

        tmp<Field<Type>> pResult(
            new Field<Type>(
                cells.size(),
                pTraits<Type>::zero
            )
        );
        Field<Type> &result=const_cast<Field<Type>&>(pResult());

        typedef VolField<Type> volField;

        word interpolationScheme("none");
        if (
            !interpolationSchemes_.found(name)
            &&
            interpolationSchemes_.found("default")
        ) {
            interpolationScheme=word(interpolationSchemes_["default"]);
        } else if (interpolationSchemes_.found(name)) {
            interpolationScheme=word(interpolationSchemes_[name]);
        }

        if (
            this->searchInMemory()
            &&
            this->mesh().foundObject<volField>(name)
        ) {
            const volField &f=this->mesh().lookupObject<volField>(name);

            autoPtr<interpolation<Type>> inter;
            inter=interpolation<Type>::New(
                interpolationScheme,
                f
            );
            forAll(result,i) {
                result[i]=inter().interpolate(
                    positions[i],
                    cells[i]
                );
            }
        } else {
            autoPtr<volField> f=this->getOrReadField<volField>(
                name,
                true // fail if not found
            );

            autoPtr<interpolation<Type>> inter=interpolation<Type>::New(
                interpolationScheme,
                f
            );
            forAll(result,i) {
                result[i]=inter().interpolate(
                    positions[i],
                    cells[i]
                );
            }
        }

        return pResult;
   }
}


// ************************************************************************* //
