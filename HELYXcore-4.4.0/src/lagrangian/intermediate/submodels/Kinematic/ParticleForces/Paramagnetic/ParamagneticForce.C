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

#include "submodels/Kinematic/ParticleForces/Paramagnetic/ParamagneticForce.H"
#include "include/demandDrivenData.H"
#include "global/constants/electromagnetic/electromagneticConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParamagneticForce<CloudType>::ParamagneticForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    HdotGradHName_
    (
        this->coeffs().template lookupOrDefault<word>("HdotGradH", "HdotGradH")
    ),
    HdotGradHInterpPtr_(nullptr),
    magneticSusceptibility_
    (
        this->coeffs().template lookup<scalar>("magneticSusceptibility")
    )
{}


template<class CloudType>
Foam::ParamagneticForce<CloudType>::ParamagneticForce
(
    const ParamagneticForce& pf
)
:
    ParticleForce<CloudType>(pf),
    HdotGradHName_(pf.HdotGradHName_),
    HdotGradHInterpPtr_(pf.HdotGradHInterpPtr_),
    magneticSusceptibility_(pf.magneticSusceptibility_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParamagneticForce<CloudType>::~ParamagneticForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParamagneticForce<CloudType>::cacheFields(const bool store)
{
    if (store)
    {
        const volVectorField& HdotGradH =
            this->mesh().template lookupObject<volVectorField>(HdotGradHName_);

        HdotGradHInterpPtr_ = interpolation<vector>::New
        (
            this->owner().solution().interpolationSchemes(),
            HdotGradH
        ).ptr();
    }
    else
    {
        deleteDemandDrivenData(HdotGradHInterpPtr_);
    }
}


template<class CloudType>
Foam::forceSuSp Foam::ParamagneticForce<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero, 0.0);

    const interpolation<vector>& HdotGradHInterp = *HdotGradHInterpPtr_;

    value.Su()=
        mass*3.0*constant::electromagnetic::mu0.value()/p.rho()
       *magneticSusceptibility_/(magneticSusceptibility_ + 3)
       *HdotGradHInterp.interpolate(p.position(), p.currentTetIndices());

    return value;
}


// ************************************************************************* //
