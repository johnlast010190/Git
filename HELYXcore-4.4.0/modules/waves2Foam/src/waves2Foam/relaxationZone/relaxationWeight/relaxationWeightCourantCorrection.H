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
    Foam::relaxationWeights::relaxationWeightCourantCorrection

Description
    This class performs an additional weighting of the relaxation zone weights,
    which are based on the local cell Courant number.

    This is based on a presentation by Sopheak Seng, Technical University of
    Denmark.

SourceFiles
    relaxationWeightCourantCorrection.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef relaxationWeightCourantCorrection_H
#define relaxationWeightCourantCorrection_H

#include "fvMesh/fvMesh.H"
#include "cfdTools/general/include/fvCFD.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationWeights
{

/*---------------------------------------------------------------------------*\
                      Class relaxationShape Declaration
\*---------------------------------------------------------------------------*/

class relaxationWeightCourantCorrection
{
private:

    // Private data

        const fvMesh& mesh_;

        Switch courantCorrection_;


    // Private Member Functions

        scalar courantNumber
        (
            const surfaceScalarField& phi,
            const surfaceScalarField& magSf,
            const surfaceScalarField& deltaCoeffs,
            const cell& cc
        ) const;


public:

    //- Runtime type information
    TypeName("relaxationWeightCourantCorrection");


    // Constructors

        //- Construct from components
        relaxationWeightCourantCorrection
        (
            const fvMesh&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        relaxationWeightCourantCorrection
        (
            const relaxationWeightCourantCorrection&
        ) = delete;


    //- Destructor
    virtual ~relaxationWeightCourantCorrection();


    // Member Functions

        void courantCorrection
        (
            const labelList&,
            scalarField&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const relaxationWeightCourantCorrection&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationWeights
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
