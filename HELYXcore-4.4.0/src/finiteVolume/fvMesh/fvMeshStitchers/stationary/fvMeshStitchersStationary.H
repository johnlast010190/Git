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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2022-2023 Engys Ltd.

Class
    Foam::fvMeshStitchers::stationary

Description
    Mesh stitcher for stationary meshes.

SourceFiles
    fvMeshStitchersStationary.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshStitchersStationary_H
#define fvMeshStitchersStationary_H

#include "fvMesh/fvMeshStitchers/fvMeshStitcher/fvMeshStitcher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshStitchers
{

/*---------------------------------------------------------------------------*\
                 Class fvMeshStitchers::stationary Declaration
\*---------------------------------------------------------------------------*/

class stationary
:
    public fvMeshStitcher
{
    // Private Member Functions

        //- Return whether or not this stitcher supports a changing mesh
        virtual bool changing() const
        {
            return false;
        }

        //- Correct the mesh fluxes following a conform operation
        virtual void conformCorrectMeshPhi(surfaceScalarField& phi);

        //- Initialise correction of the mesh fluxes
        virtual void createNonConformalCorrectMeshPhiGeometry
        (
            SurfaceFieldBoundary<label>& polyFacesBf,
            surfaceVectorField& SfSf,
            surfaceVectorField& CfSf
        );

        //- Correct the mesh fluxes following an unconform operation
        virtual void unconformCorrectMeshPhi
        (
            const SurfaceFieldBoundary<label>& polyFacesBf,
            surfaceVectorField& SfSf,
            surfaceVectorField& CfSf,
            surfaceScalarField& phi
        );


public:

    //- Runtime type information
    TypeName("stationary");


    // Constructors

        //- Construct from fvMesh
        explicit stationary(fvMesh&);


    //- Destructor
    virtual ~stationary();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshStitchers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
