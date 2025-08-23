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
    (c) 2023-2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::ReferenceFrame

Description
    Contains the dict data that defines a reference frame

SourceFiles
    referenceFrame.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "types/coordinateSystemType.H"
#include "referenceFrames/coordinateFrame.H"
#include "postDict/postDictKeys.H"
#include "hash/hasher.H"
#include "meshes/primitiveShapes/point/point.H"

#include "vtkNew.h"
#include "vtkSmartPointer.h"
#include "vtkMatrix4x4.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class ReferenceFrameInitializer
{
    const coordinateFrame* coordinateFrame_;
    label updateTime_ = -1;
public:
    const coordinateFrame &getCoordinateFrame();
    virtual ~ReferenceFrameInitializer() = default;
protected:
    const fvMesh &defaultRegion_;
    explicit ReferenceFrameInitializer(const fvMesh &defaultRegion);
    virtual const coordinateFrame* initializeCoordinateFrame() = 0;
};

struct ReferenceFrame
{
    std::shared_ptr<ReferenceFrameInitializer> initializer_;

    explicit ReferenceFrame(const fvMesh &defaultRegion);
    ReferenceFrame(const fvMesh &defaultRegion, const word &name);
    ReferenceFrame(const fvMesh &defaultRegion, const word &name, const dictionary &dictionary);

    [[nodiscard]] static bool isReferenceFrameDict(const dictionary &dictionary, const word& name);

    [[nodiscard]] vector convertLocalVectorToGlobal(vector local) const;

    [[nodiscard]] point convertLocalPointToGlobal(point local) const;

    [[nodiscard]] vector applyLocalToGlobalTransformation(vector local) const;

    [[nodiscard]] CoordinateSystemType getType() const;

    [[nodiscard]] point getOrigin() const;

    [[nodiscard]] vector e1() const;

    [[nodiscard]] vector e2() const;

    [[nodiscard]] vector e3() const;

    [[nodiscard]] vtkSmartPointer<vtkMatrix4x4> toTransformMatrix() const;

    bool operator==(const ReferenceFrame& other) const;

private:
    [[nodiscard]] const coordinateSystem& getUpdatedCoordinateSystem() const;
};

} // End namespace

// ************************************************************************* //
