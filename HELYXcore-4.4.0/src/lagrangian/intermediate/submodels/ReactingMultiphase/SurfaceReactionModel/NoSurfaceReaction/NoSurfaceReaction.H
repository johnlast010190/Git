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
    Foam::NoSurfaceReaction

Group
    grpLagrangianIntermediateSurfaceReactionSubModels

Description
    Dummy surface reaction model for 'none'

\*---------------------------------------------------------------------------*/

#ifndef NoSurfaceReaction_H
#define NoSurfaceReaction_H

#include "submodels/ReactingMultiphase/SurfaceReactionModel/SurfaceReactionModel/SurfaceReactionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class NoSurfaceReaction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoSurfaceReaction
:
    public SurfaceReactionModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoSurfaceReaction(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoSurfaceReaction(const NoSurfaceReaction<CloudType>& srm);

        //- Construct and return a clone
        virtual autoPtr<SurfaceReactionModel<CloudType>> clone() const
        {
            return autoPtr<SurfaceReactionModel<CloudType>>
            (
                new NoSurfaceReaction<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoSurfaceReaction();


    // Member Functions

        //- Flag to indicate whether model activates devolatisation model
        virtual bool active() const;

        //- Update surface reactions
        virtual scalar calculate
        (
            const scalar dt,
            const label celli,
            const scalar d,
            const scalar T,
            const scalar Tc,
            const scalar pc,
            const scalar rhoc,
            const scalar mass,
            const scalarField& YGas,
            const scalarField& YLiquid,
            const scalarField& YSolid,
            const scalarField& YMixture,
            const scalar N,
            scalarField& dMassGas,
            scalarField& dMassLiquid,
            scalarField& dMassSolid,
            scalarField& dMassSRCarrier
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "submodels/ReactingMultiphase/SurfaceReactionModel/NoSurfaceReaction/NoSurfaceReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
