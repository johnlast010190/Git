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
    (c) 2011-2013 OpenFOAM Foundation

Description
    Make instances of a function object for all known intermediate cloud types

Contributors/Copyright:
    2016 Bernhard F.W. Gschaider <bgschaid@ice-sf.at>

\*---------------------------------------------------------------------------*/

#ifndef makeFunctionObjectsForKnownClouds_H
#define makeFunctionObjectsForKnownClouds_H

#include "include/swakCloudTypes.H"

#include "submodels/CloudFunctionObjects/CloudFunctionObject/CloudFunctionObject.H"

#include "clouds/derived/basicKinematicCloud/basicKinematicCloud.H"
#include "clouds/derived/basicThermoCloud/basicThermoCloud.H"
#include "clouds/derived/basicReactingCloud/basicReactingCloud.H"
#include "clouds/derived/basicReactingMultiphaseCloud/basicReactingMultiphaseCloud.H"

#ifdef FOAM_HAS_MPPIC_CLOUD
#include "clouds/derived/basicKinematicMPPICCloud/basicKinematicMPPICCloud.H"
#define makeCloudFunctionObjectTypeMultiMPPIC(FunctionObjectType)   \
      makeCloudFunctionObjectTypeMulti(FunctionObjectType,basicKinematicMPPICCloud);
#else
#define makeCloudFunctionObjectTypeMultiMPPIC(FunctionObjectType)
#endif

#ifdef FOAM_HAS_COLLIDING_CLOUD
#include "clouds/derived/basicKinematicCollidingCloud/basicKinematicCollidingCloud.H"
#define makeCloudFunctionObjectTypeMultiColliding(FunctionObjectType)   \
    makeCloudFunctionObjectTypeMulti(FunctionObjectType,basicKinematicCollidingCloud);
#else
#define makeCloudFunctionObjectTypeMultiColliding(FunctionObjectType)
#endif

// redefinition of the original because the typedef is not unique
#define makeCloudFunctionObjectTypeMulti(SS, CloudType)             \
                                                                        \
    typedef Foam::CloudType::kinematicCloudType KCT##CloudType;         \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<KCT##CloudType>, 0);   \
                                                                        \
    Foam::CloudFunctionObject<KCT##CloudType>::                         \
    adddictionaryConstructorToTable<Foam::SS<KCT##CloudType>>          \
    add##SS##CloudType##kinematicCloudType##ConstructorToTable_;

#define makeFunctionObjectsForKnownClouds(FunctionObjectType)       \
    makeCloudFunctionObjectTypeMulti(FunctionObjectType,basicKinematicCloud); \
    makeCloudFunctionObjectTypeMultiMPPIC(FunctionObjectType);          \
    makeCloudFunctionObjectTypeMultiColliding(FunctionObjectType);      \
    makeCloudFunctionObjectTypeMulti(FunctionObjectType,basicThermoCloud); \
    makeCloudFunctionObjectTypeMulti(FunctionObjectType,basicReactingCloud); \
    makeCloudFunctionObjectTypeMulti(FunctionObjectType,basicReactingMultiphaseCloud);

#endif
