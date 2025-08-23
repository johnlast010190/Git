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
    (c) 2017 Wikki Ltd.
    (c) 2016 Engys Ltd.

Class
    Foam::leastSquaresFaVectors

Description
    Least-squares gradient scheme vectors for the Finite Area method

SourceFiles
    leastSquaresFaVectors.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresFaVectors_H
#define leastSquaresFaVectors_H

#include "meshes/MeshObject/MeshObject.H"
#include "faMesh/faMesh.H"
#include "fields/edgeFields/edgeFieldsFwd.H"
#include "primitives/Pair/labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyTopoChangeMap;

/*---------------------------------------------------------------------------*\
                     Class leastSquaresFaVectors Declaration
\*---------------------------------------------------------------------------*/

class leastSquaresFaVectors
:
    public MeshObject<faMesh, MoveableMeshObject, leastSquaresFaVectors>
{
    // Private data

        //- Least-squares gradient vectors
        mutable edgeVectorField* pVectorsPtr_;
        mutable edgeVectorField* nVectorsPtr_;


    // Private member functions

        //- Construct Least-squares gradient vectors
        void calcLeastSquaresVectors() const;


public:

    // Declare name of the class and its debug switch
    TypeName("leastSquaresFaVectors");


    // Constructors

        //- Construct given an faMesh
        explicit leastSquaresFaVectors(const faMesh&);


    // Destructor

        virtual ~leastSquaresFaVectors();


    // Member functions

        //- Return reference to owner least square vectors
        const edgeVectorField& pVectors() const;

        //- Return reference to neighbour least square vectors
        const edgeVectorField& nVectors() const;

        //- Update after mesh motion:
        //- Delete the least square vectors when the mesh moves
        virtual bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
