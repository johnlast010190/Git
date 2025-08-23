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
    Foam::movingConeTopoFvMesh

Description
    Sample topoChangerFvMesh that moves an object in x direction
    and introduces/removes layers.

SourceFiles
    movingConeTopoFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef movingConeTopoFvMesh_H
#define movingConeTopoFvMesh_H

#include "topoChangerFvMesh/topoChangerFvMesh.H"
#include "motionSolvers/motionSolver/motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                      Class movingConeTopoFvMesh Declaration
\*---------------------------------------------------------------------------*/

class movingConeTopoFvMesh
:
    public topoChangerFvMesh
{
    // Private data

        //- Motion dictionary
        dictionary motionDict_;

        //- Motion velocity amplitude
        vector motionVelAmplitude_;

        //- Motion velocity period
        scalar motionVelPeriod_;

        //- Motion velocity period
        vector curMotionVel_;

        //- Left edge
        scalar leftEdge_;

        //- Current left obstacle position
        scalar curLeft_;

        //- Current right obstacle position
        scalar curRight_;

        //- Vertex motion mask
        scalarField motionMask_;


    // Private Member Functions

        //- Add mixer zones and modifiers
        void addZonesAndModifiers();

        //- Markup motion vertices
        tmp<scalarField> vertexMarkup
        (
            const pointField& p,
            const scalar curLeft,
            const scalar curRight
        ) const;


public:

    //- Runtime type information
    TypeName("movingConeTopoFvMesh");


    // Constructors

        //- Construct from database
        explicit movingConeTopoFvMesh(const IOobject& io);

        //- Disallow default bitwise copy construction
        movingConeTopoFvMesh(const movingConeTopoFvMesh&) = delete;


    //- Destructor
    virtual ~movingConeTopoFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const movingConeTopoFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
