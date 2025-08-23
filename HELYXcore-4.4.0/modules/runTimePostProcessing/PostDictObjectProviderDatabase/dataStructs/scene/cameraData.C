/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : Dev
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
    (c) 2015-2019 OpenCFD Ltd.
    (c) 2024-2025 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::CameraData

Description
    Information about a camera from the scene
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#include "cameraData.H"

#include "types/coordinateType.H"
#include "postDict/postDictKeys.H"
#include "db/dictionary/dictionary.H"
#include "storage/referenceFrames.H"

#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

CameraData::CameraData()
{
    this->focalPoint = point(0,0,0);
    this->position = point(0,1,0);
    this->parallelProjection = false;
    this->up = vector(0,0,1);
    this->parallelScale = 0;
    this->name = "";
    this->referenceFrame_ = nullptr;
}

CameraData::CameraData(const ReferenceFrames& rfs, const dictionary& sceneCameraDict) : CameraData()
{
    readDict(rfs, sceneCameraDict);
}

void CameraData::readDict(const ReferenceFrames& rfs, const dictionary& cameraDict)
{
    focalPoint = cameraDict.lookup<point>(cameraKeys::FOCAL_POINT_KEY);
    position = cameraDict.lookup<point>(cameraKeys::POSITION_KEY);
    up = cameraDict.lookup<Vector<scalar>>(cameraKeys::UP_KEY);
    position = readPositionToCartesian(cameraDict);
    parallelProjection = cameraDict.lookupOrDefault(cameraKeys::PARALLEL_PROJECTION_KEY, false);
    if (parallelProjection)
    {
        parallelScale = cameraDict.lookup<scalar>(cameraKeys::PARALLEL_SCALE_KEY);
    }
    else
    {
        parallelScale = 1;
    }
    std::string tempName = cameraDict.lookupOrDefault<string>(cameraKeys::NAME_KEY, "");
    std::replace(tempName.begin(), tempName.end(), ' ', '_');
    name = tempName;

    word referenceFrameName = cameraDict.lookupOrDefault<word>(
        cameraKeys::REFERENCE_FRAME_KEY,
        meshObjectsKeys::GLOBAL_COORDINATE_SYSTEM_KEY
    );
    this->referenceFrame_ = &rfs.getReferenceFrame(referenceFrameName);
}

point CameraData::readPositionToCartesian(const dictionary& dict) const
{
    auto coordType = dict.lookupOrDefault(cameraKeys::COORDINATE_TYPE_KEY, CoordinateType());
    auto readPosition = dict.lookup<point>(cameraKeys::POSITION_KEY);

    switch (coordType.getValue())
    {
        case CoordinateType::Value::CARTESIAN:
        default:
            return readPosition;
        case CoordinateType::Value::SPHERICAL_XY:
            return sphericalXYToCartesian(readPosition);
        case CoordinateType::Value::SPHERICAL_XZ:
            return sphericalXZToCartesian(readPosition);
        case CoordinateType::Value::SPHERICAL_YZ:
            return sphericalYZToCartesian(readPosition);
    }
}

point CameraData::sphericalXYToCartesian(const point& sphericalXY) const
{
    point local = sphericalToLocalCartesian(sphericalXY);
    //             Azimuth = X     Cross = Y       Up = Z
    return point(local.x(), local.y(), local.z()) + focalPoint;
}

point CameraData::sphericalXZToCartesian(const point& sphericalXZ) const
{
    point local = sphericalToLocalCartesian(sphericalXZ);
    //               Cross = X       Up = Y      Azimuth = Z
    return point(local.y(), local.z(), local.x()) + focalPoint;
}

point CameraData::sphericalYZToCartesian(const point& sphericalYZ) const
{
    point local = sphericalToLocalCartesian(sphericalYZ);
    //                Up = X      Azimuth = Y     Cross = Z
    return point(local.z(), local.x(), local.y()) + focalPoint;
}

point CameraData::sphericalToLocalCartesian(const point &spherical)
{
    scalar r = spherical.x();
    scalar radTheta = M_PI * spherical.y() / 180.0;
    scalar radPhi = M_PI * spherical.z() / 180.0;

    scalar cosTheta = std::cos(radTheta);
    scalar sinTheta = std::sin(radTheta);
    scalar cosPhi   = std::cos(radPhi);
    scalar sinPhi   = std::sin(radPhi);

    // Cartesian coordinates
    scalar x = r * cosTheta * cosPhi;
    scalar y = r * sinTheta * cosPhi;
    scalar z = r * sinPhi;
    return {x, y, z};
}

CameraData CameraData::toGlobalCamera() const
{
    CameraData globalCamera;
    if (referenceFrame_)
    {
        globalCamera.focalPoint = referenceFrame_->convertLocalPointToGlobal(focalPoint);
        globalCamera.position = referenceFrame_->convertLocalPointToGlobal(position);
        globalCamera.up = referenceFrame_->convertLocalVectorToGlobal(up);
    }
    else
    {
        globalCamera.focalPoint = focalPoint;
        globalCamera.position = position;
        globalCamera.up = up;
    }
    globalCamera.parallelProjection = parallelProjection;
    globalCamera.parallelScale = parallelScale;
    globalCamera.name = name;
    globalCamera.referenceFrame_ = nullptr;

    return globalCamera;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
