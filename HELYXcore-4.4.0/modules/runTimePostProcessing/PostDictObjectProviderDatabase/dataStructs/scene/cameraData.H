/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
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

#pragma once

#include "meshes/primitiveShapes/point/point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{
class ReferenceFrames;
class ReferenceFrame;

struct CameraData
{
    point focalPoint;
    point position;
    bool parallelProjection;
    Vector<scalar> up;
    scalar parallelScale;
    word name;
    const ReferenceFrame* referenceFrame_;

    CameraData();

    CameraData(const ReferenceFrames& rfs, const dictionary& sceneCameraDict);

    void readDict(const ReferenceFrames& rfs, const dictionary& cameraDict);

    CameraData toGlobalCamera() const;

private:
    [[nodiscard]] point readPositionToCartesian(const dictionary& dict) const;

    [[nodiscard]] point sphericalXYToCartesian(const point& sphericalXY) const;
    [[nodiscard]] point sphericalXZToCartesian(const point& sphericalXZ) const;
    [[nodiscard]] point sphericalYZToCartesian(const point& sphericalYZ) const;

    static point sphericalToLocalCartesian(const point& spherical) ;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
