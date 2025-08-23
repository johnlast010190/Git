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

#include "coalCloudScatteringFactorPluginFunction.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "include/swakCloudTypes.H"

#ifdef FOAM_REACTINGCLOUD_TEMPLATED
#include "CoalCloud.H"
#else
#include "coalCloud/coalCloud.H"
#endif

namespace Foam {

defineTypeNameAndDebug(coalCloudScatteringFactorPluginFunction,0);
addNamedToRunTimeSelectionTable(FieldValuePluginFunction,coalCloudScatteringFactorPluginFunction , name, coalCloudScatteringFactor);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coalCloudScatteringFactorPluginFunction::coalCloudScatteringFactorPluginFunction(
    const FieldValueExpressionDriver &parentDriver,
    const word &name
):
    lcsScatteringFactorPluginFunction(
        parentDriver,
        name
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

autoPtr<volScalarField> coalCloudScatteringFactorPluginFunction::internalEvaluate()
{
    // pick up the first fitting class
#ifdef FOAM_REACTINGCLOUD_TEMPLATED
    tryCall(volScalarField,constThermoCoalCloud,reactingMultiphaseCloud,sigmap());
    tryCall(volScalarField,thermoCoalCloud,reactingMultiphaseCloud,sigmap());
    tryCall(volScalarField,icoPoly8ThermoCoalCloud,reactingMultiphaseCloud,sigmap());
#else
    tryCall(volScalarField,coalCloud,reactingMultiphaseCloud,sigmap());
#endif

    return lcsScatteringFactorPluginFunction::internalEvaluate();
}

// * * * * * * * * * * * * * * * Concrete implementations * * * * * * * * * //


} // namespace

// ************************************************************************* //
