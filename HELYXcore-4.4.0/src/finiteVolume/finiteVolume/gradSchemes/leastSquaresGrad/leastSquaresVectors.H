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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2010-2024 Engys Ltd.

Class
    Foam::leastSquaresVectors

Description
    Least-squares gradient scheme vectors

SourceFiles
    leastSquaresVectors.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresVectors_H
#define leastSquaresVectors_H

#include "meshes/MeshObject/MeshObject.H"
#include "fvMesh/fvMesh.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "primitives/Pair/labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class leastSquaresVectors Declaration
\*---------------------------------------------------------------------------*/

class leastSquaresVectors
:
    public MeshObject<fvMesh, MoveableMeshObject, leastSquaresVectors>
{
    // Private data

        //- Least-squares gradient vectors
        mutable surfaceVectorField* pVectorsPtr_;
        mutable surfaceVectorField* nVectorsPtr_;
        mutable volScalarField* wLayerCellsPtr_;

    // Private Member Functions

        //- Construct Least-squares gradient vectors
        void calcLeastSquaresVectors() const;

        //- Calculate cell layer indentifier
        void calcWallLayerCells(const label) const;


public:

    // Declare name of the class and its debug switch
    TypeName("leastSquaresVectors");


    // Constructors

        //- Construct given an fvMesh
        explicit leastSquaresVectors(const fvMesh&);


    //- Destructor
    virtual ~leastSquaresVectors();


    // Member functions

        //- Return reference to owner least square vectors
        const surfaceVectorField& pVectors() const;

        //- Return reference to neighbour least square vectors
        const surfaceVectorField& nVectors() const;

        //- Return reference of near wall cell marker
        const volScalarField& wallLayerCells(const label) const;

        //- Update after mesh motion:
        //  Delete the least square vectors when the mesh moves
        virtual bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
