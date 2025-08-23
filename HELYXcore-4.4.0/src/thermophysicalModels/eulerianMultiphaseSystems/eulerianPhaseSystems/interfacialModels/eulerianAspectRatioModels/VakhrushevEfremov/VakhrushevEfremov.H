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
    (c) 2014-2015 OpenFOAM Foundation

Class
    Foam::VakhrushevEfremov

Description
    Aspect ratio model of Vakhrushev and Efremov.

    Reference:
    \verbatim
        "Interpolation formula for computing the velocities of single gas
         bubbles in liquids"
        Vakhrushev, I.A. and Efremov, G.I.,
        Chemistry and Technology of Fuels and Oils
        Volume 6, Issue 5, May 1970, pp. 376-379,
    \endverbatim

SourceFiles
    VakhrushevEfremov.C

\*---------------------------------------------------------------------------*/

#ifndef VakhrushevEfremov_H
#define VakhrushevEfremov_H

#include "../eulerianAspectRatioModel/eulerianAspectRatioModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianAspectRatioModels
{

/*---------------------------------------------------------------------------*\
                       Class VakhrushevEfremov Declaration
\*---------------------------------------------------------------------------*/

class VakhrushevEfremov
:
    public eulerianAspectRatioModel
{
public:

    //- Runtime type information
    TypeName("VakhrushevEfremov");


    // Constructors

        //- Construct from a dictionary and an ordered phase pair
        VakhrushevEfremov
        (
            const dictionary& dict,
            const eulerianPhasePair& pair
        );


    //- Destructor
    virtual ~VakhrushevEfremov();


    // Member Functions

        //- Aspect ratio
        virtual tmp<volScalarField> E() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace eulerianAspectRatioModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
