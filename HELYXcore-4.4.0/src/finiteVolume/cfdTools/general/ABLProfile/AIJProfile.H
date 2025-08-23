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
    (c) 2017 Engys Ltd.

Class
    Foam::AIJProfile

Description
    AIJProfile for ABL inlet specification

SourceFiles
    AIJProfile.C

\*---------------------------------------------------------------------------*/

#ifndef AIJProfile_H
#define AIJProfile_H

#include "cfdTools/general/ABLProfile/ABLProfile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ABLProfiles
{

/*---------------------------------------------------------------------------*\
                          Class AIJProfile Declaration
\*---------------------------------------------------------------------------*/

class AIJProfile
:
    public ABLProfile
{
    // Private Data

        //- Cmu constant
        const scalar Cmu_;

        //- Reference velocity
        const vector Uref_;

        //- Reference height
        const scalar Href_;

        //- Power law exponent determined by terrain category
        const scalar alpha_;

        //- boundary layer height determined by terrain category
        const scalar zG_;

        //- tke calculated as (UI)^2 based on AIJ paper
        const bool AIJPaper_;

        //- scaling for L to compute Lxx, Lyy, Lzz from
        //  L = k*sqrt(k)/epsilon * pow(Cmu_,0.75)
        const vector fLxyz_;


public:

    //- Runtime type information
    TypeName("AIJ");

    // Constructors

        //- Null constructor
        AIJProfile(const fvPatch&);

        //- Construct from components
        AIJProfile(const fvPatch&, const dictionary&);

        //- Copy constructor
        AIJProfile(const AIJProfile&);

        //- return a clone
        virtual ABLProfile* clone() const
        {
            return new AIJProfile(*this);
        }


    //- Destructor
    virtual ~AIJProfile(){}


    // Member Functions

        //- Set boundary fields
        virtual void setProfiles
        (
            vectorField& U,
            scalarField& L,
            symmTensorField& R
        );

        //- Set boundary fields
        virtual void setProfiles
        (
            vectorField& U,
            vectorField& L,
            symmTensorField& R
        );

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const AIJProfile&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ABLProfiles
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
