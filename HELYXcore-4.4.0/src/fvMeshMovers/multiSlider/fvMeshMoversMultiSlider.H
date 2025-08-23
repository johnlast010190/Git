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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2017-2025 Engys Ltd.

Class
    Foam::fvMeshMovers:multiSlider

Description
    fvMeshMover for slider motion.

SourceFiles
    fvMeshMoversMultiSlider.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshMoversMultiSlider_H
#define fvMeshMoversMultiSlider_H

#include "fvMesh/fvMeshMovers/fvMeshMover/fvMeshMover.H"
#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshMovers
{

/*---------------------------------------------------------------------------*\
                  Class fvMeshMovers::multiSlider Declaration
\*---------------------------------------------------------------------------*/

class multiSlider
:
    public fvMeshMover
{
    // Private Data

        //- Dictionary of motion control parameters
        const dictionary coeffs_;

        //- Calculated mesh displacement
        vectorField displacement_;

        //- Specified cellZones
        labelList zoneIDs_;

        //- Points to move per cellZone
        labelListList pointIDs_;

        //- Moving patches per cellZone
        labelListList movingIDs_;

        //- Moving patches velocity or displacement per cellZone
        PtrList<Function1<vector>> motionProfile_;

        //- Whether velocity or displacement profile
        boolList velocityProfile_;

        //- Initial displacement per cellZone
        vectorField initialDisplacement_;

        //- Buffer region per cellZone
        scalarField bufferSize_;

        //- Blend type per cellZone
        wordList blendType_;


public:

    //- Runtime type information
    TypeName("multiSlider");


    // Constructors

        //- Construct from fvMesh
        multiSlider(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        multiSlider(const multiSlider&) = delete;


    //- Destructor
    ~multiSlider();


    // Member Functions

        //- Generate indirect primitive patch from moving patch ID's
        autoPtr<indirectPrimitivePatch> makePatch(const labelList& patchIDs);

        // Calculate sliding displacement by which to move the mesh
        void slidePoints(const vectorField& zoneDisplacements);

        //- Update the mesh for mesh motion
        virtual bool update() override;

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&) override
        {
            NotImplemented;
        }

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&) override
        {}

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&) override
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiSlider&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
