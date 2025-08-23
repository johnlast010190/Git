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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "global/constants/mathematical/mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::DSMCParcel<ParcelType>::constantProperties::constantProperties()
:
    mass_(0),
    d_(0)
{}


template<class ParcelType>
inline Foam::DSMCParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& dict
)
:
    mass_(dict.lookup<scalar>("mass")),
    d_(dict.lookup<scalar>("diameter")),
    internalDegreesOfFreedom_
    (
        readInt(dict.lookup("internalDegreesOfFreedom"))
    ),
    omega_(dict.lookup<scalar>("omega"))
{}


template<class ParcelType>
inline Foam::DSMCParcel<ParcelType>::DSMCParcel
(
    const polyMesh& mesh,
    const vector& position,
    const vector& U,
    const scalar Ei,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label typeId
)
:
    ParcelType(mesh, position, celli, tetFacei, tetPti),
    U_(U),
    Ei_(Ei),
    typeId_(typeId)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::DSMCParcel<ParcelType>::constantProperties::mass() const
{
    return mass_;
}


template<class ParcelType>
inline Foam::scalar Foam::DSMCParcel<ParcelType>::constantProperties::d() const
{
    return d_;
}


template<class ParcelType>
inline Foam::scalar
Foam::DSMCParcel<ParcelType>::constantProperties::sigmaT() const
{
    return constant::mathematical::pi*d_*d_;
}


template<class ParcelType>
inline Foam::direction
Foam::DSMCParcel<ParcelType>::constantProperties::internalDegreesOfFreedom()
const
{
    return internalDegreesOfFreedom_;
}


template<class ParcelType>
inline Foam::scalar
Foam::DSMCParcel<ParcelType>::constantProperties::omega() const
{
    return omega_;
}


// * * * * * * * * * * DSMCParcel Member Functions  * * * * * * * * * * //

template<class ParcelType>
inline Foam::label Foam::DSMCParcel<ParcelType>::typeId() const
{
    return typeId_;
}


template<class ParcelType>
inline const Foam::vector& Foam::DSMCParcel<ParcelType>::U() const
{
    return U_;
}


template<class ParcelType>
inline Foam::scalar Foam::DSMCParcel<ParcelType>::Ei() const
{
    return Ei_;
}


template<class ParcelType>
inline Foam::vector& Foam::DSMCParcel<ParcelType>::U()
{
    return U_;
}


template<class ParcelType>
inline Foam::scalar& Foam::DSMCParcel<ParcelType>::Ei()
{
    return Ei_;
}


// ************************************************************************* //
