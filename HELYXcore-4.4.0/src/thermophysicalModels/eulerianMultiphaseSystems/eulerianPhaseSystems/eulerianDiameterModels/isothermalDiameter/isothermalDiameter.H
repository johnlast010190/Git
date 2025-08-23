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
    Foam::eulerianDiameterModels::isothermal

Description
    Isothermal dispersed-phase particle diameter model.

SourceFiles
    isothermal.C

\*---------------------------------------------------------------------------*/

#ifndef isothermal_H
#define isothermal_H

#include "../eulerianDiameterModel/eulerianDiameterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianDiameterModels
{

/*---------------------------------------------------------------------------*\
                           Class isothermal Declaration
\*---------------------------------------------------------------------------*/

class isothermal
:
    public eulerianDiameterModel
{
    // Private data

        //- Reference diameter for the isothermal expansion
        dimensionedScalar d0_;

        //- Reference pressure for the isothermal expansion
        dimensionedScalar p0_;


public:

    //- Runtime type information
    TypeName("isothermal");


    // Constructors

        //- Construct from components
        isothermal
        (
            const dictionary& diameterProperties,
            const eulerianPhaseModel& phase
        );


    //- Destructor
    virtual ~isothermal();


    // Member Functions

        //- Return the diameter field
        virtual tmp<volScalarField> d() const;

        //- Read phaseProperties dictionary
        virtual bool read(const dictionary& phaseProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace eulerianDiameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
