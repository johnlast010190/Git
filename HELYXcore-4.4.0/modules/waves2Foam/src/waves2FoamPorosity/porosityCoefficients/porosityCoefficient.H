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
    (c) held by original author
    (c) 2024 Engys Ltd.

Class
    Foam::porosityCoefficient

Description
    Abstract base class for the computation of resistance coefficients in the
    case of porous simulations.

SourceFiles
    porosityCoefficient.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef porosityCoefficient_H
#define porosityCoefficient_H

#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class porosityCoefficient Declaration
\*---------------------------------------------------------------------------*/

class porosityCoefficient
{
protected:

    // Protected data

        const dictionary & poroProperties_;

        dimensionedVector linearCoefficient_;
        dimensionedVector quadraticCoefficient_;
        dimensionedVector KCQuadraticCoefficient_;
        dimensionedScalar scaledKC_;


    // Protected member functions

        scalar readResistancePorosity(const dictionary&) const;


public:

    //- Runtime type information
    TypeName("porosityCoefficient");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            porosityCoefficient,
            porosityCoefficient,
            (
                const dictionary & poroProp
            ),
            (poroProp)
        );


    // Constructors

        //- Construct from components
        porosityCoefficient(const dictionary &);

        //- Disallow default bitwise copy construction
        porosityCoefficient(const porosityCoefficient&) = delete;


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<porosityCoefficient> New(const dictionary &);


    //- Destructor
    virtual ~porosityCoefficient();


    // Member Functions

        // Access
        const dimensionedVector& linearCoefficient() const
        {
            return linearCoefficient_;
        };

        const dimensionedVector& quadraticCoefficient() const
        {
            return quadraticCoefficient_;
        };

        const dimensionedVector& KCquadraticCoefficinet() const
        {
        	return KCQuadraticCoefficient_;
        }

        const dimensionedScalar& scaledKC() const
        {
        	return scaledKC_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const porosityCoefficient&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
