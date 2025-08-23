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
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

template<class ParcelType>
template<class CloudType>
inline Foam::KinematicParcel<ParcelType>::TrackingData<CloudType>::TrackingData
(
    CloudType& cloud,
    trackPart part
)
:
    ParcelType::template TrackingData<CloudType>(cloud),
    rhoInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.rho()
        )
    ),
    UInterp_
    (
        interpolation<vector>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.U()
        )
    ),
    muInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.mu()
        )
    ),
    g_(cloud.g().value()),
    part_(part)
{}


template<class ParcelType>
template<class CloudType>
inline const Foam::interpolation<Foam::scalar>&
Foam::KinematicParcel<ParcelType>::TrackingData<CloudType>::rhoInterp() const
{
    return rhoInterp_();
}


template<class ParcelType>
template<class CloudType>
inline const Foam::interpolation<Foam::vector>&
Foam::KinematicParcel<ParcelType>::TrackingData<CloudType>::UInterp() const
{
    return UInterp_();
}


template<class ParcelType>
template<class CloudType>
inline const Foam::interpolation<Foam::scalar>&
Foam::KinematicParcel<ParcelType>::TrackingData<CloudType>::muInterp() const
{
    return muInterp_();
}


template<class ParcelType>
template<class CloudType>
inline const Foam::vector&
Foam::KinematicParcel<ParcelType>::TrackingData<CloudType>::g() const
{
    return g_;
}


template<class ParcelType>
template<class CloudType>
inline typename Foam::KinematicParcel<ParcelType>::template
TrackingData<CloudType>::trackPart
Foam::KinematicParcel<ParcelType>::TrackingData<CloudType>::part() const
{
    return part_;
}


template<class ParcelType>
template<class CloudType>
inline typename Foam::KinematicParcel<ParcelType>::template
TrackingData<CloudType>::trackPart&
Foam::KinematicParcel<ParcelType>::TrackingData<CloudType>::part()
{
    return part_;
}


// ************************************************************************* //
