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
    (c) 2024-2025 Engys Ltd.

Class
    Foam::fvMeshGIBChanger

Description
    Abstract base class for fvMesh GIB changers.

SourceFiles
    fvMeshGIBChanger.C
    fvMeshGIBChangerNew.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshGIBChanger_H
#define fvMeshGIBChanger_H

#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class fvMeshGIBChanger Declaration
\*---------------------------------------------------------------------------*/

class fvMeshGIBChanger
{
    // Private Data

        //- Non-const fvMesh reference to allow update
        fvMesh& mesh_;

        //- The dynamicMeshDict
        dictionary dynamicMeshDict_;


public:

    //- Runtime type information
    TypeName("fvMeshGIBChanger");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fvMeshGIBChanger,
            fvMesh,
            (fvMesh& mesh, const dictionary& dict),
            (mesh, dict)
        );


    // Constructors

        //- Construct from fvMesh
        explicit fvMeshGIBChanger(fvMesh&);

        //- Disallow default bitwise copy construction
        fvMeshGIBChanger(const fvMeshGIBChanger&) = delete;


    // Selectors

        //- Select, construct and return the fvMeshGIBChanger.
        //  If the constant/dynamicMeshDict dictionary does not exist,
        //  a dummy fvMeshGIBChangers::none object is returned.
        static autoPtr<fvMeshGIBChanger> New(fvMesh&);


    //- Destructor
    virtual ~fvMeshGIBChanger();


    // Member Functions

        //- Access the fvMesh
        fvMesh& mesh()
        {
            return mesh_;
        }

        //- Return the fvMesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return the dynamicMeshDict/GIBChanger sub-dict
        const dictionary& dict() const
        {
            return dynamicMeshDict_.subDict("GIBChanger");
        }

        //- Is mesh dynamic, i.e. might it change?
        //  Defaults to true, set to false in the fvMeshGIBChangers::none
        virtual bool dynamic() const
        {
            return true;
        }

        //- Update the mesh for changes in GIB patches
        virtual bool update() = 0;

        //- Write the changer state
        virtual bool write(const bool write = true) const
        {
            return true;
        }

        //- Return old points for mesh motion
        virtual const pointField& oldPoints() const
        {
            return mesh().polyMesh::oldPoints();
        }

        //- Correct velocity for GIB
        virtual tmp<vectorField> velocityCorrect(const vectorField& fc) const
        {
            return tmp<vectorField>(new vectorField(fc.size(), vector::zero));
        }

        //- faceZone initialization using the stl information
        virtual void updateInit(const word&)
        {};

        //- Write GIB geometry
        virtual void writeGeometry(const fileName& surfaceName)
        {};


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvMeshGIBChanger&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
