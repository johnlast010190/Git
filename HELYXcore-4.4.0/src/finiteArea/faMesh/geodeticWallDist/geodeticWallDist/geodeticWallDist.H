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
    (c) 2016-2024 Engys Ltd.
    (c) 2015 OpenFOAM Foundation

Class
    Foam::geodeticWallDist

Description
    Interface to run-time selectable methods to calculate the distance-to-wall
    and normal-to-wall fields.

    Example of the wallDist specification in fvSchemes:
    \verbatim
        wallDist
        {
            method Poisson;

            // Optional entry enabling the calculation
            // of the normal-to-wall field
            nRequired false;
        }
    \endverbatim

SeeAlso
    Foam::faPatchDistMethod::Poisson
    Foam::faPatchDistMethod::advectionDiffusion

SourceFiles
    geodeticWallDist.C

\*---------------------------------------------------------------------------*/

#ifndef geodeticWallDist_H
#define geodeticWallDist_H

#include "meshes/MeshObject/MeshObject.H"
#include "faMesh/geodeticWallDist/faPatchDistMethods/faPatchDistMethod/faPatchDistMethod.H"
#include "fields/areaFields/areaFields.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class geodeticWallDist Declaration
\*---------------------------------------------------------------------------*/

class geodeticWallDist
:
    public MeshObject<faMesh, UpdateableMeshObject, geodeticWallDist>
{
    // Private data

        //- Set of patch IDs
        const labelHashSet patchIDs_;

        //- island faces
        mutable labelList islandFaces_;

        //- Name for the patch set, e.g. "wall"
        const word patchTypeName_;

        //- Run-time selected method to generate the distance-to-wall field
        mutable autoPtr<faPatchDistMethod> pdm_;

        //- Distance-to-wall field
        mutable areaScalarField y_;

        //- Flag to indicate if the distance-to-wall field is required
        mutable bool nRequired_;

        //- Normal-to-wall field
        mutable tmp<areaVectorField> n_;

    // Private Member Functions

        //- Read PDM method (set to Poisson if not present)
        const dictionary& setPDM
        (
            const faMesh& mesh,
            const word& patchTypeName
        );

        //- Construct the normal-to-wall field as required
        void constructn() const;


public:

    // Declare name of the class and its debug switch
    ClassName("geodeticWallDist");


    // Constructors

        //- Construct from mesh and optional patch type name
        geodeticWallDist
        (
            const faMesh& mesh,
            const labelHashSet&,
            const word& patchTypeName = "wall"
        );

        //- Construct from mesh, patch IDs and optional patch type name
        geodeticWallDist
        (
            const faMesh& mesh,
            const labelHashSet& patchIDs,
            const labelList&,
            const word& patchTypeName = "patch"
        );

        //- Disallow default bitwise copy construction
        geodeticWallDist(const geodeticWallDist&) = delete;


    //- Destructor
    virtual ~geodeticWallDist();

    // Member Functions

        //- Return the patchIDs
        const labelHashSet& patchIDs() const
        {
            return patchIDs_;
        }

        //- Return reference to cached distance-to-wall field
        const areaScalarField& y() const
        {
            return y_;
        }

        //- Return reference to cached normal-to-wall field
        const areaVectorField& n() const;

        //- Update the y-field when the mesh moves
        virtual bool movePoints();

        //- Update the y-field when the mesh changes
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&)
        {
            NotImplemented;
        }

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&)
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const geodeticWallDist&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
