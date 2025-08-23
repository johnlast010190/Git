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

Class
    Foam::relaxationWeights::relaxationWeightThirdOrderPolynomial

Description
    An object used to create a weighting function in the relaxation zone, which
    has a specific polynomial behaviour. See eq. 6.18 in

    A. P. Engsig-Karup
    Unstructured Nodal DG-FEM Solution of High-Order Boussinesq-Type Equations
    PhD Thesis
    Technical University of Denmark

SourceFiles
    relaxationWeightThirdOrderPolynomial.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef relaxationWeightThirdOrderPolynomial_H
#define relaxationWeightThirdOrderPolynomial_H

#include "relaxationZone/relaxationWeight/relaxationWeight.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationWeights
{

/*---------------------------------------------------------------------------*\
                           Class relaxationShapeRectangular Declaration
\*---------------------------------------------------------------------------*/

class relaxationWeightThirdOrderPolynomial
:
    public relaxationWeight
{

private:

    // Private data

protected:

    // Protected Member functions
        virtual void computeWeights
        (
            const labelList&,
            const scalarField&,
            scalarField&
        );

public:

    //- Runtime type information
    TypeName("relaxationWeightThirdOrderPolynomial");

    // Constructors

        //- from components
        relaxationWeightThirdOrderPolynomial
        (
            const word& ,
            const fvMesh& mesh_
        );


    // Destructor

        ~relaxationWeightThirdOrderPolynomial()
        {}


    // Member Functions

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationWeights
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
