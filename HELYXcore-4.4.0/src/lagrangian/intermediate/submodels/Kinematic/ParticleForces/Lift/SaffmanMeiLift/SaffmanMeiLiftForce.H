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
    (c) 2012-2016 OpenFOAM Foundation

Class
    Foam::SaffmanMeiLiftForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Saffman-Mei particle lift force model applicable to spherical particles.

SourceFiles
    SaffmanMeiLiftForce.C

\*---------------------------------------------------------------------------*/

#ifndef SaffmanMeiLiftForce_H
#define SaffmanMeiLiftForce_H

#include "submodels/Kinematic/ParticleForces/Lift/LiftForce/LiftForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class SaffmanMeiLiftForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SaffmanMeiLiftForce
:
    public LiftForce<CloudType>
{
protected:

    // Protected Member Functions

        //- Calculate the lift coefficient
        virtual scalar Cl
        (
            const typename CloudType::parcelType& p,
            const vector& curlUc,
            const scalar Re,
            const scalar muc
        ) const;


public:

    //- Runtime type information
    TypeName("SaffmanMeiLiftForce");


    // Constructors

        //- Construct from mesh
        SaffmanMeiLiftForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& forceType = typeName
        );

        //- Construct copy
        SaffmanMeiLiftForce(const SaffmanMeiLiftForce& lf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new SaffmanMeiLiftForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~SaffmanMeiLiftForce();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "submodels/Kinematic/ParticleForces/Lift/SaffmanMeiLift/SaffmanMeiLiftForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
