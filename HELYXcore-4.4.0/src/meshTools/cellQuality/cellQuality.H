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
    (c) 2011 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::cellQuality

Description
    Class calculates cell quality measures.

SourceFiles
    cellQuality.C

\*---------------------------------------------------------------------------*/

#ifndef cellQuality_H
#define cellQuality_H

#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class cellQuality Declaration
\*---------------------------------------------------------------------------*/

class cellQuality
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;


public:

    // Constructors

        //- Construct from mesh
        cellQuality(const polyMesh& mesh);

        //- Disallow default bitwise copy construction
        cellQuality(const cellQuality&) = delete;


    // Member Functions

        //- Return cell non-orthogonality
        tmp<scalarField> nonOrthogonality() const;

        //- Return cell skewness
        tmp<scalarField> skewness() const;

        //- Return face non-orthogonality
        tmp<scalarField> faceNonOrthogonality() const;

        //- Return face skewness
        tmp<scalarField> faceSkewness() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cellQuality&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
