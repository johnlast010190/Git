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

#include "lcsMassFractionPluginFunction.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "include/swakCloudTypes.H"

#include "clouds/derived/basicKinematicCloud/basicKinematicCloud.H"
#ifdef FOAM_REACTINGCLOUD_TEMPLATED
#include "clouds/derived/basicThermoCloud/basicThermoCloud.H"
#include "BasicReactingCloud.H"
#include "BasicReactingMultiphaseCloud.H"
#else
#include "clouds/derived/basicReactingCloud/basicReactingCloud.H"
#include "clouds/derived/basicReactingMultiphaseCloud/basicReactingMultiphaseCloud.H"
#endif

namespace Foam {

defineTypeNameAndDebug(lcsMassFractionPluginFunction,0);
addNamedToRunTimeSelectionTable(FieldValuePluginFunction,lcsMassFractionPluginFunction , name, lcsMassFraction);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lcsMassFractionPluginFunction::lcsMassFractionPluginFunction(
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

autoPtr<volScalarField> lcsMassFractionPluginFunction::internalEvaluate()
{
    // pick up the first fitting class
    tryCall(volScalarField,basicKinematicCloud,kinematicCloud,alpha());
#ifdef FOAM_REACTINGCLOUD_TEMPLATED
    tryCall(volScalarField,basicThermoCloud,thermoCloud,alpha());
    tryCall(volScalarField,constThermoReactingCloud,reactingCloud,alpha());
    tryCall(volScalarField,thermoReactingCloud,reactingCloud,alpha());
    tryCall(volScalarField,icoPoly8ThermoReactingCloud,reactingCloud,alpha());
    tryCall(volScalarField,constThermoReactingMultiphaseCloud,reactingMultiphaseCloud,alpha());
    tryCall(volScalarField,thermoReactingMultiphaseCloud,reactingMultiphaseCloud,alpha());
    tryCall(volScalarField,icoPoly8ThermoReactingMultiphaseCloud,reactingMultiphaseCloud,alpha());
#else
    tryCall(volScalarField,swakFluidThermoCloudType,thermoCloud,alpha());
    tryCall(volScalarField,basicReactingCloud,reactingCloud,alpha());
    tryCall(volScalarField,basicReactingMultiphaseCloud,reactingMultiphaseCloud,alpha());
#endif

    return autoPtr<volScalarField>();
}

void lcsMassFractionPluginFunction::doEvaluation()
{
    autoPtr<volScalarField> palpha=internalEvaluate();

    noCloudFound(palpha);

    result().setObjectResult(palpha);
}

// * * * * * * * * * * * * * * * Concrete implementations * * * * * * * * * //


} // namespace

// ************************************************************************* //
