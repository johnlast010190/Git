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
    2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "include/swak.H"

#ifdef  FOAM_FV_HAS_SMOOTH_SWEEP_SPREAD

#include "fvcSweepPluginFunction.H"
#include "FieldValueExpressionDriver.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "finiteVolume/fvc/fvcSmooth/fvcSmooth.H"

namespace Foam {

defineTypeNameAndDebug(fvcSweepPluginFunction,1);
addNamedToRunTimeSelectionTable(FieldValuePluginFunction, fvcSweepPluginFunction , name, fvcSweep);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fvcSweepPluginFunction::fvcSweepPluginFunction(
    const FieldValueExpressionDriver &parentDriver,
    const word &name
):
    FieldValuePluginFunction(
        parentDriver,
        name,
        word("volScalarField"),
        string(
            "originalField internalField volScalarField"
            ",alphaField internalField volScalarField"
            ",nLayers primitive label"
            ",alphaDiff_default=0.2 primitive scalar"
        )
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fvcSweepPluginFunction::doEvaluation()
{
    autoPtr<volScalarField> pResult(
        new volScalarField(
            IOobject(
                "sweep",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless, 0),
            "zeroGradient"
        )
    );
    volScalarField &result=pResult();

    result.forceAssign(field_());

    fvc::sweep(
        result,
        alpha_,
        nLayers_,
        alphaDiff_
    );

    this->result().setObjectResult(pResult);
}

void fvcSweepPluginFunction::setArgument(
    label index,
    const string &content,
    const CommonValueExpressionDriver &driver
) {
    assert(index==0 || index==1);

    if (index==0) {
        this->field_.set(
            new volScalarField(
                dynamic_cast<const FieldValueExpressionDriver &>(
                    driver
                ).getResult<volScalarField>()
            )
        );
    } else {
        this->alpha_.set(
            new volScalarField(
                dynamic_cast<const FieldValueExpressionDriver &>(
                    driver
                ).getResult<volScalarField>()
            )
        );
    }
}

void fvcSweepPluginFunction::setArgument(
    label index,
    const scalar &val
)
{
    assert(index==3);

    alphaDiff_=val;
}

void fvcSweepPluginFunction::setArgument(
    label index,
    const label &val
)
{
    assert(index==2);
    nLayers_=val;
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace

#endif

// ************************************************************************* //
