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

#include "lcsRhoEffPluginFunction.H"

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

defineTypeNameAndDebug(lcsRhoEffPluginFunction,0);
addNamedToRunTimeSelectionTable(FieldValuePluginFunction,lcsRhoEffPluginFunction , name, lcsRhoEff);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lcsRhoEffPluginFunction::lcsRhoEffPluginFunction(
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

autoPtr<volScalarField> lcsRhoEffPluginFunction::internalEvaluate()
{
    // pick up the first fitting class
    tryCall(volScalarField,basicKinematicCloud,kinematicCloud,rhoEff());
#ifdef FOAM_REACTINGCLOUD_TEMPLATED
    tryCall(volScalarField,basicThermoCloud,thermoCloud,rhoEff());
    tryCall(volScalarField,constThermoReactingCloud,reactingCloud,rhoEff());
    tryCall(volScalarField,thermoReactingCloud,reactingCloud,rhoEff());
    tryCall(volScalarField,icoPoly8ThermoReactingCloud,reactingCloud,rhoEff());
    tryCall(volScalarField,constThermoReactingMultiphaseCloud,reactingMultiphaseCloud,rhoEff());
    tryCall(volScalarField,thermoReactingMultiphaseCloud,reactingMultiphaseCloud,rhoEff());
    tryCall(volScalarField,icoPoly8ThermoReactingMultiphaseCloud,reactingMultiphaseCloud,rhoEff());
#else
    tryCall(volScalarField,swakFluidThermoCloudType,thermoCloud,rhoEff());
    tryCall(volScalarField,basicReactingCloud,reactingCloud,rhoEff());
    tryCall(volScalarField,basicReactingMultiphaseCloud,reactingMultiphaseCloud,rhoEff());
#endif

    return autoPtr<volScalarField>();
}

void lcsRhoEffPluginFunction::doEvaluation()
{

    autoPtr<volScalarField> prhoEff=internalEvaluate();

    noCloudFound(prhoEff);

    result().setObjectResult(prhoEff);
}

// * * * * * * * * * * * * * * * Concrete implementations * * * * * * * * * //


} // namespace

// ************************************************************************* //
