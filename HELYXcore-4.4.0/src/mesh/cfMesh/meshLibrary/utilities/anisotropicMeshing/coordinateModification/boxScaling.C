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
    (c) Creative Fields, Ltd.
    (c) 2024 Engys Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

\*---------------------------------------------------------------------------*/

#include "utilities/anisotropicMeshing/coordinateModification/boxScaling.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "meshes/boundBox/boundBox.H"
#include "meshes/primitiveShapes/plane/plane.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(boxScaling, 0);
addToRunTimeSelectionTable(coordinateModification, boxScaling, dictionary);


// * * * * * * * * * * * * * * Private member functions* * * * * * * * * * * //

void boxScaling::calculateBndBox()
{
    pMin_ = centre_ - 0.5 * lengthVec_;
    pMax_ = centre_ + 0.5 * lengthVec_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

boxScaling::boxScaling()
:
    coordinateModification(),
    centre_(),
    lengthVec_(0.0, 0.0, 0.0),
    scaleVec_(1.0, 1.0, 1.0),
    pMin_(),
    pMax_()
{
    calculateBndBox();
}


boxScaling::boxScaling
(
    const word& name,
    const point& centre,
    const scalar lengthX,
    const scalar lengthY,
    const scalar lengthZ,
    const scalar scaleX,
    const scalar scaleY,
    const scalar scaleZ
)
:
    coordinateModification(),
    centre_(centre),
    lengthVec_(lengthX, lengthY, lengthZ),
    scaleVec_(scaleX, scaleY, scaleZ),
    pMin_(),
    pMax_()
{
    calculateBndBox();
    setName(name);
}

boxScaling::boxScaling
(
    const word& name,
    const dictionary& dict
)
:
    coordinateModification(name, dict)
{
    this->operator=(dict);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

point boxScaling::origin() const
{
    return centre_;
}


void boxScaling::translateAndModifyObject(const vector& disp)
{
    centre_ += disp;

    for (direction i=0;i<vector::nComponents;++i)
        lengthVec_[i] /= scaleVec_[i];

    calculateBndBox();
}


vector boxScaling::displacement(const point& p) const
{
    vector disp;

    for (direction i=0;i<vector::nComponents;++i)
    {
        const scalar dispVec = lengthVec_[i] * ((1.0/scaleVec_[i]) - 1.0);
        const scalar t = ((p[i] - pMin_[i]) / lengthVec_[i]);

        const scalar tBnd = Foam::max(0.0, Foam::min(t, 1.0));

        disp[i] = tBnd * dispVec;
    }

    return disp;
}


vector boxScaling::backwardDisplacement(const point& p) const
{
    vector disp;

    for (direction i=0;i<vector::nComponents;++i)
    {
        const scalar dispVec = lengthVec_[i] * (scaleVec_[i] - 1.0);

        const scalar t = ((p[i] - pMin_[i]) / lengthVec_[i]);

        const scalar tBnd = Foam::max(0.0, Foam::min(t, 1.0));

        disp[i] = tBnd * dispVec;
    }

    return disp;
}


bool boxScaling::combiningPossible() const
{
    return true;
}


void boxScaling::boundingPlanes(PtrList<plane>&pl) const
{
    pl.setSize(6);
    label counter(0);
    if (Foam::mag(scaleVec_.x() - 1.0) > VSMALL)
    {
        pl.set(counter++, new plane(pMin_, vector(1, 0, 0)));
        pl.set(counter++, new plane(pMax_, vector(1, 0, 1)));
    }

    if (Foam::mag(scaleVec_.y() - 1.0) > VSMALL)
    {
        pl.set(counter++, new plane(pMin_, vector(0, 1, 0)));
        pl.set(counter++, new plane(pMax_, vector(0, 1, 0)));
    }

    if (Foam::mag(scaleVec_.z() - 1.0) > VSMALL)
    {
        pl.set(counter++, new plane(pMin_, vector(0, 0, 1)));
        pl.set(counter++, new plane(pMax_, vector(0, 0, 1)));
    }

    pl.setSize(counter);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dictionary boxScaling::dict(bool /*ignoreType*/) const
{
    dictionary dict;

    dict.add("type", type());

    dict.add("centre", centre_);
    dict.add("lengthX", lengthVec_.x());
    dict.add("lengthY", lengthVec_.y());
    dict.add("lengthZ", lengthVec_.z());

    dict.add("scaleX", scaleVec_.x());
    dict.add("scaleY", scaleVec_.y());
    dict.add("scaleZ", scaleVec_.z());

    return dict;
}


void boxScaling::write(Ostream& os) const
{
    os  << " type:   " << type()
        << " centre: " << centre_
        << " lengthX: " << lengthVec_.x()
        << " lengthY: " << lengthVec_.y()
        << " lengthZ: " << lengthVec_.z()
        << " scaleX:  " << scaleVec_.x()
        << " scaleY:  " << scaleVec_.y()
        << " scaleZ:  " << scaleVec_.z()
        << endl;
}


void boxScaling::writeDict(Ostream& os, bool subDict) const
{
    if (subDict)
    {
        os << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    // only write type for derived types
    if (type() != typeName_())
    {
        os.writeEntry("type", type());
    }

    os.writeEntry("centre", centre_);
    os.writeEntry("lengthX", lengthVec_.x());
    os.writeEntry("lengthY", lengthVec_.y());
    os.writeEntry("lengthZ", lengthVec_.z());
    os.writeEntry("scaleX", scaleVec_.x());
    os.writeEntry("scaleY", scaleVec_.y());
    os.writeEntry("scaleZ", scaleVec_.z());

    if (subDict)
    {
        os.endBlock();
    }
}


void boxScaling::operator=(const dictionary& d)
{
    // allow as embedded sub-dictionary "coordinateSystem"
    const dictionary& dict = d.found(typeName_()) ? d.subDict(typeName_()) : d;

    // Unspecified centre is (0 0 0)
    centre_ = dict.lookup<point>("centre");

    // Specify lengthX
    lengthVec_.x() = dict.lookup<scalar>("lengthX");

    // Specify lengthY
    lengthVec_.y() = dict.lookup<scalar>("lengthY");

    // Specify lengthZ
    lengthVec_.z() = dict.lookup<scalar>("lengthZ");

    // Specify scaleX
    scaleVec_.x() = dict.lookupOrDefault<scalar>("scaleX", 1.0);

    // Specify scaleY
    scaleVec_.y() = dict.lookupOrDefault<scalar>("scaleY", 1.0);

    // Specify scaleX
    scaleVec_.z() = dict.lookupOrDefault<scalar>("scaleZ", 1.0);

    calculateBndBox();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Ostream& boxScaling::operator<<(Ostream& os) const
{
    os << "name " << name() << nl;
    write(os);
    return os;
}

Ostream& operator<<(Ostream& os, const boxScaling& bs)
{
    return bs.operator<<(os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
