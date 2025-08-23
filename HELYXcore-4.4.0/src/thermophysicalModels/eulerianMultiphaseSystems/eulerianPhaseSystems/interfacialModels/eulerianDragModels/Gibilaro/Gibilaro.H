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
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::eulerianDragModels::Gibilaro

Description
    H, Enwald, E. Peirano, A-E Almstedt
    'Eulerian Two-Phase Flow Theory Applied to Fluidization'
    Int. J. Multiphase Flow, Vol. 22, Suppl, pp. 21-66 (1996)
    Eq. 106, p. 43

SourceFiles
    Gibilaro.C

\*---------------------------------------------------------------------------*/

#ifndef Gibilaro_H
#define Gibilaro_H

#include "../eulerianDragModel/eulerianDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class eulerianPhasePair;

namespace eulerianDragModels
{

/*---------------------------------------------------------------------------*\
                           Class Gibilaro Declaration
\*---------------------------------------------------------------------------*/

class Gibilaro
:
    public eulerianDragModel
{
public:

    //- Runtime type information
    TypeName("Gibilaro");


    // Constructors

        //- Construct from a dictionary and a phase pair
        Gibilaro
        (
            const dictionary& dict,
            const eulerianPhasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~Gibilaro();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace eulerianDragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
