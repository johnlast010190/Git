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
    (c) 2013-2016 OpenFOAM Foundation

Class
    Foam::ParticleStressModels::exponential

Group
    grpLagrangianIntermediateMPPICParticleStressSubModels

Description
    Exponential inter-particle stress model of the same form as used in
    twoPhaseEulerFoam

SourceFiles
    exponential.C

\*---------------------------------------------------------------------------*/

#ifndef exponential_H
#define exponential_H

#include "submodels/MPPIC/ParticleStressModels/ParticleStressModel/ParticleStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ParticleStressModels
{

/*---------------------------------------------------------------------------*\
                        Class exponential Declaration
\*---------------------------------------------------------------------------*/

class exponential
:
    public ParticleStressModel
{
    // Private data

        //- Pre-exponential factor
        scalar preExp_;

        //- Maximum limit of the exponential
        scalar expMax_;

        //- Front coefficient
        scalar g0_;


public:

    //- Runtime type information
    TypeName("exponential");


    //- Constructors

        //- Construct from components
        exponential(const dictionary& dict);

        //- Construct copy
        exponential(const exponential& hc);

        //- Clone
        virtual autoPtr<ParticleStressModel> clone() const
        {
            return autoPtr<ParticleStressModel>
            (
                new exponential(*this)
            );
        }


    //- Destructor
    virtual ~exponential();


    //- Member Functions

        //- Collision stress
        tmp<Field<scalar>> tau
        (
            const Field<scalar>& alpha,
            const Field<scalar>& rho,
            const Field<scalar>& uRms
        ) const;

        //- Collision stress derivaive w.r.t. the volume fraction
        tmp<Field<scalar>> dTaudTheta
        (
            const Field<scalar>& alpha,
            const Field<scalar>& rho,
            const Field<scalar>& uRms
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ParticleStressModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
