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

#include "lcsMomentumSourcePluginFunction.H"

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

defineTypeNameAndDebug(lcsMomentumSourcePluginFunction,0);
addNamedToRunTimeSelectionTable(FieldValuePluginFunction,lcsMomentumSourcePluginFunction , name, lcsMomentumSource);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lcsMomentumSourcePluginFunction::lcsMomentumSourcePluginFunction(
    const FieldValueExpressionDriver &parentDriver,
    const word &name
):
    LagrangianCloudSourcePluginFunction(
        parentDriver,
        name,
        "volVectorField"
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

autoPtr<lcsMomentumSourcePluginFunction::dimVectorField>
lcsMomentumSourcePluginFunction::internalEvaluate()
{
    // pick up the first fitting class
#ifdef FOAM_KINEMATICCLOUD_OLD_STYLE
    tryCall(dimVectorField,basicKinematicCloud,kinematicCloud,SU());
    tryCall(dimVectorField,basicThermoCloud,thermoCloud,SU());
    tryCall(dimVectorField,constThermoReactingCloud,reactingCloud,SU());
    tryCall(dimVectorField,thermoReactingCloud,reactingCloud,SU());
    tryCall(dimVectorField,icoPoly8ThermoReactingCloud,reactingCloud,SU());
    tryCall(dimVectorField,constThermoReactingMultiphaseCloud,reactingMultiphaseCloud,SU());
    tryCall(dimVectorField,thermoReactingMultiphaseCloud,reactingMultiphaseCloud,SU());
    tryCall(dimVectorField,icoPoly8ThermoReactingMultiphaseCloud,reactingMultiphaseCloud,SU());
#else
    tryCall(dimVectorField,basicKinematicCloud,kinematicCloud,UTrans());
    tryCall(dimVectorField,swakFluidThermoCloudType,thermoCloud,UTrans());
    tryCall(dimVectorField,basicReactingCloud,reactingCloud,UTrans());
    tryCall(dimVectorField,basicReactingMultiphaseCloud,reactingMultiphaseCloud,UTrans());
#endif

    return autoPtr<dimVectorField>();
}

void lcsMomentumSourcePluginFunction::doEvaluation()
{
    autoPtr<dimVectorField> pSU=internalEvaluate();

    noCloudFound(pSU);

    const dimVectorField &SU=pSU();

    autoPtr<volVectorField> pSource(
        new volVectorField(
            IOobject(
                cloudName()+"MomentumSource",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
#ifdef FOAM_KINEMATICCLOUD_OLD_STYLE
            SU.dimensions(),
#else
            SU.dimensions()/(dimTime*dimVolume),
#endif
            "zeroGradient"
        )
    );

#ifdef FOAM_KINEMATICCLOUD_OLD_STYLE
    pSource->internalField()=SU.field();
#else
#ifdef FOAM_NO_DIMENSIONEDINTERNAL_IN_GEOMETRIC
    const_cast<vectorField&>(pSource->internalField().field())
#else
    pSource->internalField()
#endif
        =SU.field()/(mesh().time().deltaT().value()*mesh().V());
#endif

    result().setObjectResult(pSource);
}

// * * * * * * * * * * * * * * * Concrete implementations * * * * * * * * * //


} // namespace

// ************************************************************************* //
