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
    (c) 2020-2024 Engys Ltd.

Class
    Foam::Function1Types::AIJWindKineticEnergyProfile

Description
    AIJWindKineticEnergyProfile for ABL inlet specification
    Usage:
    \verbatim
    <entryName> <AIJWindKineticEnergyProfileFunction>;
    <entryName>Coeffs
    {
        Cmu 0.09;
        Uref 5;
        alpha 0.2;
        Href 10;
        zG 10;
    }
    \endverbatim

    Turbulent intensity, k, is based on the estimation equation for
    vertical profile of intensity recommended by AIJ

    \f[
    i = \frac{\sigma(z)}{UMean} = 0.1*\frac{z}{z_G}^(-alpha-0.05)
    \f]
    \f[
    k = 3/2*sqr(i*|UMean|)
    \f]

    where z_G is the turbulent boundary layer height determined by terrain
    type.

    The additional consistentTKE flag allows computation of tke as (UI)^2 instead
    of 3.2*(UI)^2, according to AIJ Paper

See also
    Foam::Function1

SourceFiles
    AIJWindKineticEnergyProfile.C

\*---------------------------------------------------------------------------*/

#ifndef AIJWindKineticEnergyProfile_H
#define AIJWindKineticEnergyProfile_H

#include "primitives/functions/Function1/Function1/Function1.H"
#include "primitives/functions/Function1/AIJWindProfile/AIJWindProfile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                Class AIJWindKineticEnergyProfile Declaration
\*---------------------------------------------------------------------------*/

class AIJWindKineticEnergyProfile
:
    public AIJWindProfile
{
private:

    //- Turbulent kinetic energy calculated according to AIJ as (IU)^2
    bool AIJPaper_;


public:

    // Runtime type information
    TypeName("AIJWindKineticEnergyProfile");


    // Constructors

        //- Construct from entry name and dictionary
        AIJWindKineticEnergyProfile
        (
            const word& entryName,
            const dictionary& dict
        );

        //- Copy constructor
        AIJWindKineticEnergyProfile(const AIJWindKineticEnergyProfile& wp);

        //- Construct and return a clone
        virtual tmp<Function1<scalar>> clone() const override
        {
            return
                tmp<Function1<scalar>>
                (
                    new AIJWindKineticEnergyProfile(*this)
                );
        }


       //- Destructor
       virtual ~AIJWindKineticEnergyProfile();


    // Member Functions

        //- Return value for coordinate z
        virtual scalar value(const scalar z) const override;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const AIJWindKineticEnergyProfile&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
