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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "referenceFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

ReferenceFrameInitializer::ReferenceFrameInitializer(const fvMesh &defaultRegion)
:
    coordinateFrame_(nullptr),
    defaultRegion_(defaultRegion)
{};

const coordinateFrame &ReferenceFrameInitializer::getCoordinateFrame()
{
    label currentTime = defaultRegion_.time().timeIndex();
    if (updateTime_ != currentTime)
    {
        this->coordinateFrame_ = initializeCoordinateFrame();
        this->coordinateFrame_->updateState();
        this->updateTime_ = currentTime;
    }
    return *this->coordinateFrame_;
}

class GlobalReferenceFrameInitializer : public ReferenceFrameInitializer
{


public:
    explicit GlobalReferenceFrameInitializer(const fvMesh &defaultRegion) : ReferenceFrameInitializer(defaultRegion)
    {
        getCoordinateFrame();
    };

    const coordinateFrame *initializeCoordinateFrame() override
    {
        return coordinateFrame::globalFrame(defaultRegion_);
    }

    ~GlobalReferenceFrameInitializer() override = default;
};

class NamedReferenceFrameInitializer : public ReferenceFrameInitializer
{
    const word name_;

public:
    NamedReferenceFrameInitializer(const fvMesh &defaultRegion, const word &name)
        :
        ReferenceFrameInitializer(defaultRegion),
        name_(name)
    {
        getCoordinateFrame();
    };

    const coordinateFrame *initializeCoordinateFrame() override
    {
        return &coordinateFrame::New(defaultRegion_, name_);
    }

    ~NamedReferenceFrameInitializer() override = default;
};

class DictionaryReferenceFrameInitializer : public ReferenceFrameInitializer
{
    const word name_;
    const dictionary &dictionary_;

public:
    DictionaryReferenceFrameInitializer(const fvMesh &defaultRegion, const word &name, const dictionary &dictionary)
        :
        ReferenceFrameInitializer(defaultRegion),
        name_(name),
        dictionary_(dictionary)
    {
        getCoordinateFrame();
    };

    const coordinateFrame *initializeCoordinateFrame() override
    {
        return coordinateFrame::New(defaultRegion_, name_, dictionary_);
    }

    ~DictionaryReferenceFrameInitializer() override = default;
};

ReferenceFrame::ReferenceFrame(const fvMesh &defaultRegion)
    :
    initializer_(new GlobalReferenceFrameInitializer(defaultRegion))
{};

ReferenceFrame::ReferenceFrame(const fvMesh &defaultRegion, const word &name)
    :
    initializer_(new NamedReferenceFrameInitializer(defaultRegion, name))
{};

ReferenceFrame::ReferenceFrame(const fvMesh &defaultRegion, const word &name, const dictionary &dictionary)
    :
    initializer_(new DictionaryReferenceFrameInitializer(defaultRegion, name, dictionary))
{};

bool ReferenceFrame::isReferenceFrameDict(const dictionary &dictionary, const word& name)
{
    if (!dictionary.isDict(name)) return false;

    word type = dictionary.subDict(name).lookupOrDefault<>("type", word("unknown"));
    auto& referenceFrameTypes = coordinateFrame::dictionaryConstructorTable_();
    return referenceFrameTypes.find(type) != referenceFrameTypes.end();
}

vector ReferenceFrame::convertLocalVectorToGlobal(vector local) const
{
    return applyLocalToGlobalTransformation(local);
}

point ReferenceFrame::convertLocalPointToGlobal(point local) const
{
    point origin = getUpdatedCoordinateSystem().origin();
    return applyLocalToGlobalTransformation(local) + origin;
}

vector ReferenceFrame::applyLocalToGlobalTransformation(vector local) const
{
    const coordinateSystem &coorSys = getUpdatedCoordinateSystem();
    return coorSys.e1() * local.x() + coorSys.e2() * local.y() + coorSys.e3() * local.z();
}

CoordinateSystemType ReferenceFrame::getType() const
{
    return CoordinateSystemType(getUpdatedCoordinateSystem().type());
}

point ReferenceFrame::getOrigin() const
{
    return getUpdatedCoordinateSystem().origin();
}

vector ReferenceFrame::e1() const
{
    return getUpdatedCoordinateSystem().e1();
}

vector ReferenceFrame::e2() const
{
    return getUpdatedCoordinateSystem().e2();
}

vector ReferenceFrame::e3() const
{
    return getUpdatedCoordinateSystem().e3();
}

vtkSmartPointer<vtkMatrix4x4> ReferenceFrame::toTransformMatrix() const
{
    vtkNew<vtkMatrix4x4> matrix;
    matrix->Identity();

    const coordinateSystem &coorSys = getUpdatedCoordinateSystem();

    vector e1 = coorSys.e1();
    vector e2 = coorSys.e2();
    vector e3 = coorSys.e3();
    point origin = coorSys.origin();
    for (int row = 0; row < 3; row++)
    {
        matrix->SetElement(row, 0, e1[row]);
        matrix->SetElement(row, 1, e2[row]);
        matrix->SetElement(row, 2, e3[row]);
        matrix->SetElement(row, 3, origin[row]);
    }
    return matrix;
}

bool ReferenceFrame::operator==(const ReferenceFrame &other) const
{
    return &this->initializer_->getCoordinateFrame() == &other.initializer_->getCoordinateFrame();
}

const coordinateSystem &ReferenceFrame::getUpdatedCoordinateSystem() const
{
    return initializer_->getCoordinateFrame().coorSys();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
