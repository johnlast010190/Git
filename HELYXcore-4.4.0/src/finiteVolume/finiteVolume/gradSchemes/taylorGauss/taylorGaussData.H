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
    (c) 2016-2024 Engys Ltd.

Class
    Foam::taylorGaussData

Description
    TaylorGauss gradient scheme - compute inverse 'pseudo volume' tensors
    derived from Taylor expansion.

SourceFiles
    taylorGaussData.C

\*---------------------------------------------------------------------------*/

#ifndef taylorGaussData_H
#define taylorGaussData_H

#include "meshes/MeshObject/MeshObject.H"
#include "fvMesh/fvMesh.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "primitives/Pair/labelPair.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class taylorGaussData Declaration
\*---------------------------------------------------------------------------*/

class taylorGaussData
:
    public MeshObject<fvMesh, MoveableMeshObject, taylorGaussData>
{
    // Private data

        //- geometric inverse pseudo-volume tensor
        mutable tensorField* invPseudoVolPtr_;

        //- field specific inverse pseudo-volume tensors
        mutable HashPtrTable<Field<tensor>, word, string::hash>*
            fieldInvPseudoVolPtrs_;

        //- Pseudo-Volume condition number
        //  Currently it is taken as max/min eigenvalue
        //  Close to 1 is healthy - large can make the system ill-conditioned
        mutable scalarField* cAPtr_;

    // Private Member Functions

        //- constrain
        void constrainSolution(tensorField&) const;

        //- Construct pseudo-volume tensors
        tmp<tensorField> pseudoVolumes() const;

        //- Calculate field specific pseudo-volume modifications
        template<class Type>
        tmp<tensorField> pseudoVolumes
        (
            VolField<Type>
        ) const;

        //- Calculate condition number
        tmp<scalarField> calcConditionNumber() const;

        //- Clear data
        void clearData();

public:

    // Declare name of the class and its debug switch
    TypeName("taylorGaussData");

    //- Cut off GG-TG for cells with bad conditioning for pseudoVol inversion
        static float condNuThreshold_;

    // Constructors

        //- Construct given an fvMesh
        explicit taylorGaussData(const fvMesh&);


    //- Destructor
    virtual ~taylorGaussData();


    // Member functions

        //- Return reference to inverse geometric pseudo volumes
        const tensorField& invPseudoVolumes() const;

        //- Return reference to condition number
        const scalarField& conditionNumber() const;

        //- Return reference to inverse field plus geometric pseudo volumes
        const tensorField& invPseudoVolumes(const volScalarField&) const;

        //- Return reference to inverse  field plus geometric pseudo volumes
        const tensorField& invPseudoVolumes(const volVectorField&) const;

        //- Delete the pseudo volume vectors when the mesh moves
        virtual bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
