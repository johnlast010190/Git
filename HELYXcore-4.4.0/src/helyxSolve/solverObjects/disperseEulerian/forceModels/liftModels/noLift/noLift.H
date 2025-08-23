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
    (c) 2014 OpenFOAM Foundation

Class
    Foam::decoupledEulerian::noLift

Description
    Constant coefficient lift model.

SourceFiles
    noLift.C

\*---------------------------------------------------------------------------*/

#ifndef noLift_H
#define noLift_H

#include "solverObjects/disperseEulerian/forceModels/liftModels/liftModel/liftModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phase;

namespace decoupledEulerian
{

/*---------------------------------------------------------------------------*\
                           Class noLift Declaration
\*---------------------------------------------------------------------------*/

class noLift
:
    public liftModel
{
    // Private data


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from a dictionary and a phase pair
        noLift
        (
            const dictionary& dict,
            const phase& phase,
            const bool registerObject
        );


    //- Destructor
    virtual ~noLift();


    // Member Functions

        //- Lift coefficient
        virtual tmp<volScalarField> Cl() const;

        //- Lift force
        virtual tmp<volVectorField> F() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace decoupledEulerian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
