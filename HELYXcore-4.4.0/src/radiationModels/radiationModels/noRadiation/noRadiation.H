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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::radiationModels::noRadiation

Description
    No radiation - does nothing to energy equation source terms
    (returns zeros)

SourceFiles
    noRadiation.C

\*---------------------------------------------------------------------------*/

#ifndef noRadiation_H
#define noRadiation_H

#include "radiationModels/radiationModel/radiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{

/*---------------------------------------------------------------------------*\
                         Class noRadiation Declaration
\*---------------------------------------------------------------------------*/

class noRadiation
:
    public radiationModel
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noRadiation(const volScalarField& T, dimensionedScalar TRef);

        //- Construct from components
        noRadiation
        (
            const dictionary& dict,
            const volScalarField& T,
            dimensionedScalar TRef
        );

        //- Disallow default bitwise copy construction
        noRadiation(const noRadiation&) = delete;


    //- Destructor
    virtual ~noRadiation();


    // Member functions

        //- Main update/correction routine
        virtual bool correct() override;

        //- Solve radiation equation(s)
        virtual void calculate() override;

        //- Read radiationProperties dictionary
        virtual bool read() override;

        //- Source term component (for power of T^4)
        virtual tmp<volScalarField> Rp() const override;

        //- Source term component (constant)
        virtual tmp<volScalarField::Internal> Ru() const override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const noRadiation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
