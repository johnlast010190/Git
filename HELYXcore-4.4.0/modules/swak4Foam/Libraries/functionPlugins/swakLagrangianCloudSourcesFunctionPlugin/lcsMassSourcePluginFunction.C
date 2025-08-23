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

#include "lcsMassSourcePluginFunction.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "include/swakCloudTypes.H"

#ifdef FOAM_REACTINGCLOUD_TEMPLATED
#include "BasicReactingCloud.H"
#include "BasicReactingMultiphaseCloud.H"
#else
#include "clouds/derived/basicReactingCloud/basicReactingCloud.H"
#include "clouds/derived/basicReactingMultiphaseCloud/basicReactingMultiphaseCloud.H"
#endif

namespace Foam {

defineTypeNameAndDebug(lcsMassSourcePluginFunction,0);
addNamedToRunTimeSelectionTable(FieldValuePluginFunction,lcsMassSourcePluginFunction , name, lcsMassSource);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lcsMassSourcePluginFunction::lcsMassSourcePluginFunction(
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

autoPtr<lcsMassSourcePluginFunction::dimScalarField> lcsMassSourcePluginFunction::internalEvaluate()
{
    // pick up the first fitting class
#ifdef FOAM_REACTINGCLOUD_TEMPLATED
    tryCall(dimScalarField,constThermoReactingCloud,reactingCloud,Srho());
    tryCall(dimScalarField,thermoReactingCloud,reactingCloud,Srho());
    tryCall(dimScalarField,icoPoly8ThermoReactingCloud,reactingCloud,Srho());
    tryCall(dimScalarField,constThermoReactingMultiphaseCloud,reactingMultiphaseCloud,Srho());
    tryCall(dimScalarField,thermoReactingMultiphaseCloud,reactingMultiphaseCloud,Srho());
    tryCall(dimScalarField,icoPoly8ThermoReactingMultiphaseCloud,reactingMultiphaseCloud,Srho());
#else
    tryCall(dimScalarField,basicReactingCloud,reactingCloud,Srho());
    tryCall(dimScalarField,basicReactingMultiphaseCloud,reactingMultiphaseCloud,Srho());
#endif

    return autoPtr<dimScalarField>();
}

void lcsMassSourcePluginFunction::doEvaluation()
{

    autoPtr<dimScalarField> pSrho=internalEvaluate();

    noCloudFound(pSrho);

    const dimScalarField &Srho=pSrho();

    autoPtr<volScalarField> pSource(
        new volScalarField(
            IOobject(
                cloudName()+"MassSource",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            Srho.dimensions(),
            "zeroGradient"
        )
    );

#ifdef FOAM_NO_DIMENSIONEDINTERNAL_IN_GEOMETRIC
    const_cast<scalarField&>(pSource->internalField().field())
#else
    pSource->internalField()
#endif
        =Srho.field();

    result().setObjectResult(pSource);
}

// * * * * * * * * * * * * * * * Concrete implementations * * * * * * * * * //


} // namespace

// ************************************************************************* //
