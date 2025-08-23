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

Class
    Foam::ParamagneticForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Calculates particle paramagnetic (magnetic field) force

SourceFiles
    ParamagneticForceI.H
    ParamagneticForce.C

\*---------------------------------------------------------------------------*/

#ifndef ParamagneticForce_H
#define ParamagneticForce_H

#include "submodels/Kinematic/ParticleForces/ParticleForce/ParticleForce.H"
#include "interpolation/interpolation/interpolation/interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                       Class ParamagneticForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParamagneticForce
:
    public ParticleForce<CloudType>
{
    // Private data

        //- Name of paramagnetic field strength field - default = "HdotGradH"
        const word HdotGradHName_;

        //- HdotGradH interpolator
        const interpolation<vector>* HdotGradHInterpPtr_;

        //- Magnetic susceptibility of particle
        const scalar magneticSusceptibility_;


public:

    //- Runtime type information
    TypeName("paramagnetic");


    // Constructors

        //- Construct from mesh
        ParamagneticForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        ParamagneticForce(const ParamagneticForce& gf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new ParamagneticForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParamagneticForce();


    // Member Functions

        // Access

            //- Return the name of paramagnetic field strength field
            const word& HdotGradHName() const;

            //- Return the magnetic susceptibility of particle
            scalar magneticSusceptibility() const;


        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

            //- Calculate the non-coupled force
            virtual forceSuSp calcNonCoupled
            (
                const typename CloudType::parcelType& p,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "submodels/Kinematic/ParticleForces/Paramagnetic/ParamagneticForceI.H"

#ifdef NoRepository
    #include "submodels/Kinematic/ParticleForces/Paramagnetic/ParamagneticForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
