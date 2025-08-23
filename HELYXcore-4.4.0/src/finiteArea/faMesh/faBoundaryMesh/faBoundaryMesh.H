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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Engys Ltd.

Class
    faBoundaryMesh

Description
    Finite area boundary mesh

SourceFiles
    faBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef faBoundaryMesh_H
#define faBoundaryMesh_H

#include "faMesh/faPatches/faPatch/faPatchList.H"
#include "matrices/lduMatrix/lduAddressing/lduInterface/lduInterfacePtrsList.H"
#include "primitives/strings/lists/wordList.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "db/regIOobject/regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class faMesh;

/*---------------------------------------------------------------------------*\
                      Class faBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class faBoundaryMesh
:
    public faPatchList,
    public regIOobject
{
    // private data

        //- Reference to mesh
        const faMesh& mesh_;


public:

    //- Runtime type information
    TypeName("faBoundaryMesh");


    // Constructors

        //- Construct from dictionary
        faBoundaryMesh
        (
            const IOobject& io,
            const faMesh& fam
        );

        //- Construct given size
        faBoundaryMesh
        (
            const IOobject& io,
            const faMesh& fam,
            const label size
        );

        //- Disallow default bitwise copy construction
        faBoundaryMesh(const faBoundaryMesh&) = delete;


    // Member functions

        // Access

            //- Calculate the geometry for the patches
            //  (transformation tensors etc.)
            void calcGeometry();

            //- Return the mesh reference
            const faMesh& mesh() const;

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            lduInterfacePtrsList interfaces() const;

            //- Return a list of patch types
            wordList types() const;

            //- Return a list of patch names
            wordList names() const;

            //- Find patch index given a name
            label findPatchID(const word& patchName) const;

            //- Return patch index for a given edge label
            label whichPatch(const label edgeIndex) const;

            //- Find patch indices given a name
            labelList findIndices(const keyType&, const bool useGroups) const;

            //- Check boundary definition
            bool checkDefinition(const bool report = false) const;

        // Edit

            //- Correct faBoundaryMesh after moving points
            void movePoints(const pointField&);

            //- Correct faBoundaryMesh after topology update
            void topoChange();

            //- writeData member function required by regIOobject
            bool writeData(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const faBoundaryMesh&) = delete;


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const faBoundaryMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
