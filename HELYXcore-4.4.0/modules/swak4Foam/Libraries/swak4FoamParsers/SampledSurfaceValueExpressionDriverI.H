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
    2010-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceMesh/surfaceMesh.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchField.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchField.H"
#include "interpolations/primitivePatchInterpolation/primitivePatchInterpolation.H"

#include <cassert>

namespace Foam {
    template<class Type,class volField,class surfField>
    tmp<Field<Type>>
    SampledSurfaceValueExpressionDriver::sampleOrInterpolateInternal(
        const word &name,
        bool oldTime
    )
    {
        if (
            this->hasVariable(name)
            &&
            this->variable(name).valueType()==pTraits<Type>::typeName
        ) {
            return this->variable(name).getResult<Type>(true);
        } else if (
            this->isGlobal<Type>(name,false)
        ) {
            const ExpressionResult &globalVar=this->lookupGlobal(name);
            if (!globalVar.isSingleValue()) {
                return const_cast<ExpressionResult&>(
                    globalVar
                ).getResult<Type>(true);
            }
        }

        autoPtr<Field<Type>> result;

        if (
            this->searchInMemory()
            &&
            this->mesh().foundObject<volField>(name)
        ) {
            if (interpolate_) {
                autoPtr<interpolation<Type>> inter(
                    interpolation<Type>::New(
                        interpolationType_,
                        (
                            !oldTime ?
                            this->mesh().lookupObject<volField>(name) :
                            this->mesh().lookupObject<volField>(name).oldTime()
                        )
                    )
                );

                result.set(
                    new Field<Type>(
                        this->theSurface_.interpolate(
                            inter()
                        )
                    )
                );
            } else {
                result.set(
                    new Field<Type>(
                        this->theSurface_.sample(
                            (
                                !oldTime ?
                                this->mesh().lookupObject<volField>(name) :
                                this->mesh().lookupObject<volField>(name).oldTime()
                            )
                        )
                    )
                );
            }
        } else {
            autoPtr<volField> f=this->getOrReadField<volField>(
                name,
                true, // fail if not found
                oldTime
            );

            if (interpolate_) {
                autoPtr<interpolation<Type>> inter(
                    interpolation<Type>::New(
                        interpolationType_,
                        (
                            !oldTime ?
                            f() :
                            f->oldTime()
                        )
                    )
                );

                result.set(
                    new Field<Type>(
                        this->theSurface_.interpolate(
                            inter()
                        )
                    )
                );
            } else {
                result.set(
                    new Field<Type>(
                        this->theSurface_.sample(
                            (
                                !oldTime ?
                                f() :
                                f->oldTime()
                            )
                        )
                    )
                );
            }
        }

        if (!result.valid()) {
            FatalErrorIn("SampledSurfaceValueExpressionDriver::sampleOrInterpolateInternal(const word &name)")
                << "No valid result"
                    << endl
                    << exit(FatalError);
        }

        if (debug) {
            Pout<< "SampledSurfaceValueExpressionDriver::sampleOrInterpolateInternal()"
                << "Type: " << pTraits<Type>::typeName
                << " name: " << name
                << " size: " << this->size() << " result: "
                << result().size()
                << endl;
        }

        return tmp<Field<Type>>(result.ptr());
    }

}


// ************************************************************************* //
