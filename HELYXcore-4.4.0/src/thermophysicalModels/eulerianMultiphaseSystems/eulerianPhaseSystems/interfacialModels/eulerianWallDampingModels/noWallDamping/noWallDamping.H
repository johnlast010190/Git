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
    (c) 2015 OpenFOAM Foundation

Class
    Foam::eulerianWallDampingModels::noWallDamping

Description

SourceFiles
    noWallDamping.C

\*---------------------------------------------------------------------------*/

#ifndef noWallDamping_H
#define noWallDamping_H

#include "../eulerianWallDampingModel/eulerianWallDampingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class eulerianPhasePair;

namespace eulerianWallDampingModels
{

/*---------------------------------------------------------------------------*\
                           Class noWallDamping Declaration
\*---------------------------------------------------------------------------*/

class noWallDamping
:
    public eulerianWallDampingModel
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noWallDamping
        (
            const dictionary& dict,
            const eulerianPhasePair& pair
        );


    //- Destructor
    virtual ~noWallDamping();


    // Member Functions

        //- Return damped coefficient
        virtual tmp<volScalarField> damp
        (
            const tmp<volScalarField>&
        ) const;

        //- Return damped force
        virtual tmp<volVectorField> damp
        (
            const tmp<volVectorField>&
        ) const;

        //- Return damped face force
        virtual tmp<surfaceScalarField> damp
        (
            const tmp<surfaceScalarField>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace eulerianWallDampingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
