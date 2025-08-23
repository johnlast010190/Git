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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.
    (c) 2017 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "parcels/Templates/KinematicParcel/KinematicParcel.H"
#include "db/IOstreams/IOstreams.H"
#include "db/IOobjects/IOField/IOField.H"
#include "Cloud/Cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::KinematicParcel<ParcelType>::propertyList_ =
    Foam::KinematicParcel<ParcelType>::propertyList();

template<class ParcelType>
Foam::string Foam::KinematicParcel<ParcelType>::propertyTypes_ =
    Foam::KinematicParcel<ParcelType>::propertyTypes();

template<class ParcelType>
std::size_t Foam::KinematicParcel<ParcelType>::sizeofFields() const
{
    return reinterpret_cast<const char *>(&rhoc_) - reinterpret_cast<const char *>(&active_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    active_(false),
    typeId_(0),
    nParticle_(0.0),
    d_(0.0),
    dTarget_(0.0),
    U_(Zero),
    rho_(0.0),
    age_(0.0),
    tTurb_(0.0),
    UTurb_(Zero),
    parcelID_(0),
    parcelForces_(Zero),
    rhoc_(0.0),
    Uc_(Zero),
    muc_(0.0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            active_ = readBool(is);
            typeId_ = readLabel(is);
            nParticle_ = readScalar(is);
            d_ = readScalar(is);
            dTarget_ = readScalar(is);
            is >> U_;
            rho_ = readScalar(is);
            age_ = readScalar(is);
            parcelID_ = readLabel(is);
            is >> parcelForces_;
            tTurb_ = readScalar(is);
            is >> UTurb_;
        }
        else
        {
            is.read(reinterpret_cast<char*>(&active_), sizeofFields());
        }
    }

    is.check(FUNCTION_NAME);
}


template<class ParcelType>
template<class CloudType>
void Foam::KinematicParcel<ParcelType>::readFields(CloudType& c)
{
    const bool read = c.size();

    ParcelType::readFields(c);

    IOField<label> active
    (
        c.fieldIOobject("active", IOobject::MUST_READ),
        read
    );
    c.checkFieldIOobject(c, active);

    IOField<label> typeId
    (
        c.fieldIOobject("typeId", IOobject::MUST_READ),
        read
    );
    c.checkFieldIOobject(c, typeId);

    IOField<scalar> nParticle
    (
        c.fieldIOobject("nParticle", IOobject::MUST_READ),
        read
    );
    c.checkFieldIOobject(c, nParticle);

    IOField<scalar> d
    (
        c.fieldIOobject("d", IOobject::MUST_READ),
        read
    );
    c.checkFieldIOobject(c, d);

    IOField<scalar> dTarget
    (
        c.fieldIOobject("dTarget", IOobject::MUST_READ),
        read
    );
    c.checkFieldIOobject(c, dTarget);

    IOField<vector> U
    (
        c.fieldIOobject("U", IOobject::MUST_READ),
        read
    );
    c.checkFieldIOobject(c, U);

    IOField<scalar> rho
    (
        c.fieldIOobject("rho", IOobject::MUST_READ),
        read
    );
    c.checkFieldIOobject(c, rho);

    IOField<scalar> age
    (
        c.fieldIOobject("age", IOobject::MUST_READ),
        read
    );
    c.checkFieldIOobject(c, age);

    IOField<label> parcelID
    (
        c.fieldIOobject("parcelID", IOobject::MUST_READ),
        read
    );
    c.checkFieldIOobject(c, parcelID);

    IOField<vector> parcelForces
    (
        c.fieldIOobject("parcelForces", IOobject::MUST_READ),
        read
    );
    c.checkFieldIOobject(c, parcelForces);

    IOField<scalar> tTurb
    (
        c.fieldIOobject("tTurb", IOobject::MUST_READ),
        read
    );
    c.checkFieldIOobject(c, tTurb);

    IOField<vector> UTurb
    (
        c.fieldIOobject("UTurb", IOobject::MUST_READ),
        read
    );
    c.checkFieldIOobject(c, UTurb);

    label i = 0;

    forAllIter(typename CloudType, c, iter)
    {
        KinematicParcel<ParcelType>& p = iter();

        p.active_ = active[i];
        p.typeId_ = typeId[i];
        p.nParticle_ = nParticle[i];
        p.d_ = d[i];
        p.dTarget_ = dTarget[i];
        p.U_ = U[i];
        p.rho_ = rho[i];
        p.age_ = age[i];
        p.parcelID_ = parcelID[i];
        p.parcelForces_ = parcelForces[i];
        p.tTurb_ = tTurb[i];
        p.UTurb_ = UTurb[i];

        i++;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::KinematicParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np = c.size();

    IOField<label> active(c.fieldIOobject("active", IOobject::NO_READ), np);
    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::NO_READ), np);
    IOField<scalar> nParticle
    (
        c.fieldIOobject("nParticle", IOobject::NO_READ),
        np
    );
    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<scalar> dTarget(c.fieldIOobject("dTarget", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<scalar> rho(c.fieldIOobject("rho", IOobject::NO_READ), np);
    IOField<scalar> age(c.fieldIOobject("age", IOobject::NO_READ), np);
    IOField<label> parcelID
    (
        c.fieldIOobject("parcelID", IOobject::NO_READ),
        np
    );
    IOField<vector> parcelForces(c.fieldIOobject("parcelForces", IOobject::NO_READ), np);
    IOField<scalar> tTurb(c.fieldIOobject("tTurb", IOobject::NO_READ), np);
    IOField<vector> UTurb(c.fieldIOobject("UTurb", IOobject::NO_READ), np);

    label i = 0;

    forAllConstIter(typename CloudType, c, iter)
    {
        const KinematicParcel<ParcelType>& p = iter();

        active[i] = p.active();
        typeId[i] = p.typeId();
        nParticle[i] = p.nParticle();
        d[i] = p.d();
        dTarget[i] = p.dTarget();
        U[i] = p.U();
        rho[i] = p.rho();
        age[i] = p.age();
        parcelID[i] = p.parcelID();
        parcelForces[i] = p.parcelForces();
        tTurb[i] = p.tTurb();
        UTurb[i] = p.UTurb();

        i++;
    }

    const bool write = np > 0;

    active.write(write);
    typeId.write(write);
    nParticle.write(write);
    d.write(write);
    dTarget.write(write);
    U.write(write);
    rho.write(write);
    age.write(write);
    parcelID.write(write);
    parcelForces.write(write);
    tTurb.write(write);
    UTurb.write(write);
}


template<class ParcelType>
template<class CloudType>
void Foam::KinematicParcel<ParcelType>::writeObjects
(
    const CloudType& c,
    objectRegistry& obr
)
{
    ParcelType::writeObjects(c, obr);

    label np = c.size();

    IOField<label>& active(cloud::createIOField<label>("active", np, obr));
    IOField<label>& typeId(cloud::createIOField<label>("typeId", np, obr));
    IOField<scalar>& nParticle
    (
        cloud::createIOField<scalar>("nParticle", np, obr)
    );
    IOField<scalar>& d(cloud::createIOField<scalar>("d", np, obr));
    IOField<scalar>& dTarget(cloud::createIOField<scalar>("dTarget", np, obr));
    IOField<vector>& U(cloud::createIOField<vector>("U", np, obr));
    IOField<scalar>& rho(cloud::createIOField<scalar>("rho", np, obr));
    IOField<scalar>& age(cloud::createIOField<scalar>("age", np, obr));
    IOField<scalar>& tTurb(cloud::createIOField<scalar>("tTurb", np, obr));
    IOField<vector>& UTurb(cloud::createIOField<vector>("UTurb", np, obr));

    label i = 0;

    forAllConstIter(typename CloudType, c, iter)
    {
        const KinematicParcel<ParcelType>& p = iter();

        active[i] = p.active();
        typeId[i] = p.typeId();
        nParticle[i] = p.nParticle();
        d[i] = p.d();
        dTarget[i] = p.dTarget();
        U[i] = p.U();
        rho[i] = p.rho();
        age[i] = p.age();
        tTurb[i] = p.tTurb();
        UTurb[i] = p.UTurb();

        i++;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const KinematicParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << bool(p.active())
            << token::SPACE << p.typeId()
            << token::SPACE << p.nParticle()
            << token::SPACE << p.d()
            << token::SPACE << p.dTarget()
            << token::SPACE << p.U()
            << token::SPACE << p.rho()
            << token::SPACE << p.age()
            << token::SPACE << p.parcelID()
            << token::SPACE << p.parcelForces()
            << token::SPACE << p.tTurb()
            << token::SPACE << p.UTurb();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.active_),
            p.sizeofFields()
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
