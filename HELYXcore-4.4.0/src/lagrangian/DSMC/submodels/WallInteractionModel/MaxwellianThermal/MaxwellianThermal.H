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
    Foam::MaxwellianThermal

Description
    Wall interaction setting microscopic velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature

\*---------------------------------------------------------------------------*/

#ifndef MaxwellianThermal_H
#define MaxwellianThermal_H

#include "submodels/WallInteractionModel/WallInteractionModel/WallInteractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class MaxwellianThermal Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MaxwellianThermal
:
    public WallInteractionModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("MaxwellianThermal");


    // Constructors

        //- Construct from dictionary
        MaxwellianThermal
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~MaxwellianThermal();


    // Member Functions

        //- Apply wall correction
        virtual void correct
        (
            typename CloudType::parcelType& p,
            const wallPolyPatch& wpp
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "submodels/WallInteractionModel/MaxwellianThermal/MaxwellianThermal.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
