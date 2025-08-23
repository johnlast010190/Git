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
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "lcsEnthalpySourcePluginFunction.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "include/swakCloudTypes.H"

#ifdef FOAM_REACTINGCLOUD_TEMPLATED
#include "clouds/derived/basicThermoCloud/basicThermoCloud.H"
#include "BasicReactingCloud.H"
#include "BasicReactingMultiphaseCloud.H"
#else
#include "clouds/derived/basicReactingCloud/basicReactingCloud.H"
#include "clouds/derived/basicReactingMultiphaseCloud/basicReactingMultiphaseCloud.H"
#endif

namespace Foam {

defineTypeNameAndDebug(lcsEnthalpySourcePluginFunction,0);
addNamedToRunTimeSelectionTable(FieldValuePluginFunction,lcsEnthalpySourcePluginFunction , name, lcsEnthalpySource);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lcsEnthalpySourcePluginFunction::lcsEnthalpySourcePluginFunction(
    const FieldValueExpressionDriver &parentDriver,
    const word &name
):
    LagrangianCloudSourcePluginFunction(
        parentDriver,
        name,
        "volScalarField"
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

autoPtr<lcsEnthalpySourcePluginFunction::dimScalarField>
lcsEnthalpySourcePluginFunction::internalEvaluate()
{
    // pick up the first fitting class
#ifdef FOAM_REACTINGCLOUD_TEMPLATED
    tryCall(dimScalarField,basicThermoCloud,thermoCloud,Sh());
    tryCall(dimScalarField,constThermoReactingCloud,reactingCloud,Sh());
    tryCall(dimScalarField,thermoReactingCloud,reactingCloud,Sh());
    tryCall(dimScalarField,icoPoly8ThermoReactingCloud,reactingCloud,Sh());
    tryCall(dimScalarField,constThermoReactingMultiphaseCloud,reactingMultiphaseCloud,Sh());
    tryCall(dimScalarField,thermoReactingMultiphaseCloud,reactingMultiphaseCloud,Sh());
    tryCall(dimScalarField,icoPoly8ThermoReactingMultiphaseCloud,reactingMultiphaseCloud,Sh());
#else
    tryCall(dimScalarField,swakFluidThermoCloudType,thermoCloud,hsTrans());
    tryCall(dimScalarField,basicReactingCloud,reactingCloud,hsTrans());
    tryCall(dimScalarField,basicReactingMultiphaseCloud,reactingMultiphaseCloud,hsTrans());
#endif

    return autoPtr<dimScalarField>();
}

void lcsEnthalpySourcePluginFunction::doEvaluation()
{
    autoPtr<dimScalarField> pSh=internalEvaluate();

    noCloudFound(pSh);

    const dimScalarField &Sh=pSh();

    autoPtr<volScalarField> pSource(
        new volScalarField(
            IOobject(
                cloudName()+"EnthalpySource",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
#ifdef FOAM_THERMOCLOUD_OLD_STYLE
            Sh.dimensions(),
#else
            Sh.dimensions()/(dimTime*dimVolume),
#endif
            "zeroGradient"
        )
    );

#ifdef FOAM_THERMOCLOUD_OLD_STYLE
    pSource->internalField()=Sh.field();
#else
#ifdef FOAM_NO_DIMENSIONEDINTERNAL_IN_GEOMETRIC
    const_cast<scalarField&>(pSource->internalField().field())
#else
    pSource->internalField()
#endif
        =Sh.field()/(mesh().V()*mesh().time().deltaT().value());
#endif

    result().setObjectResult(pSource);
}

// * * * * * * * * * * * * * * * Concrete implementations * * * * * * * * * //


} // namespace

// ************************************************************************* //
