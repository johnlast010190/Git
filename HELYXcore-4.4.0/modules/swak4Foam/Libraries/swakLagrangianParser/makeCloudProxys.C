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
    (c) 1991-2010 OpenCFD Ltd.

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "CloudProxyForReactingMultiphaseParcel.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "include/swakCloudTypes.H"

#ifdef FOAM_REACTINGCLOUD_TEMPLATED
#include "BasicReactingMultiphaseCloud.H"
#include "BasicReactingCloud.H"
#else
#include "clouds/derived/basicReactingMultiphaseCloud/basicReactingMultiphaseCloud.H"
#include "clouds/derived/basicReactingCloud/basicReactingCloud.H"
#endif

#include "clouds/derived/basicKinematicCloud/basicKinematicCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    addCloudProxyToTable(CloudProxyForKinematicParcel,basicKinematicCloud);

    addCloudProxyToTable(CloudProxyForThermoParcel,swakFluidThermoCloudType);

#ifndef FOAM_BASICKINEMATICCLOUD_HAS_NO_KINEMATIC_TYPE
    addKinematicCloudProxyToTable(CloudProxyForKinematicParcel,basicKinematicCloud);
    addKinematicCloudProxyToTable(CloudProxyForThermoParcel,swakFluidThermoCloudType);
#endif

#ifdef FOAM_REACTINGCLOUD_TEMPLATED
    addCloudProxyToTable(CloudProxyForReactingParcel,thermoReactingCloud);
    addCloudProxyToTable(CloudProxyForReactingParcel,constThermoReactingCloud);
    addCloudProxyToTable(CloudProxyForReactingParcel,icoPoly8ThermoReactingCloud);

    addCloudProxyToTable(CloudProxyForReactingMultiphaseParcel,thermoReactingMultiphaseCloud);
    addCloudProxyToTable(CloudProxyForReactingMultiphaseParcel,constThermoReactingMultiphaseCloud);
    addCloudProxyToTable(CloudProxyForReactingMultiphaseParcel,icoPoly8ThermoReactingMultiphaseCloud);
#else
    addCloudProxyToTable(CloudProxyForReactingParcel,basicReactingCloud);

    addCloudProxyToTable(CloudProxyForReactingMultiphaseParcel,basicReactingMultiphaseCloud);

#ifndef FOAM_BASICKINEMATICCLOUD_HAS_NO_KINEMATIC_TYPE
    addKinematicCloudProxyToTable(CloudProxyForReactingParcelNoComposition,basicReactingCloud);
    addKinematicCloudProxyToTable(CloudProxyForReactingMultiphaseParcelNoComposition,basicReactingMultiphaseCloud);
#endif

#endif


} // namespace end

// ************************************************************************* //
