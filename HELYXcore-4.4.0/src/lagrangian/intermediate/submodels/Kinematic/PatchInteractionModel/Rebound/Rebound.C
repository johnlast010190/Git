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
    (c) 2016 OpenCFD Ltd.
    (c) 2011-2012 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "submodels/Kinematic/PatchInteractionModel/Rebound/Rebound.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::Rebound<CloudType>::Rebound
(
    const dictionary& dict,
    CloudType& cloud
)
:
    PatchInteractionModel<CloudType>(dict, cloud, typeName),
    UFactor_(this->coeffDict().template lookup<scalar>("UFactor"))
{}


template<class CloudType>
Foam::Rebound<CloudType>::Rebound(const Rebound<CloudType>& pim)
:
    PatchInteractionModel<CloudType>(pim),
    UFactor_(pim.UFactor_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::Rebound<CloudType>::~Rebound()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::Rebound<CloudType>::correct
(
    typename CloudType::parcelType& p,
    const polyPatch& pp,
    bool& keepParticle,
    const scalar trackFraction,
    const tetIndices& tetIs
)
{
    vector& U = p.U();

    keepParticle = true;
    p.active(true);

    vector nw;
    vector Up;

    this->owner().patchData(p, pp, trackFraction, tetIs, nw, Up);

    // Calculate motion relative to patch velocity
    U -= Up;

    scalar Un = U & nw;

    if (Un > 0.0)
    {
        U -= UFactor_*2.0*Un*nw;
    }

    // Return velocity to global space
    U += Up;

    return true;
}


// ************************************************************************* //
