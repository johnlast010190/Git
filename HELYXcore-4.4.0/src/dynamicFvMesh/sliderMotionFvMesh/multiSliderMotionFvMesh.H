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
    (c) 2017-2024 Engys Ltd.

Class
    Foam::multiSliderMotionFvMesh

Description
    Solid-body motion of the mesh specified by a run-time selectable
    motion function.

SourceFiles
    multiSliderMotionFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef multiSliderMotionFvMesh_H
#define multiSliderMotionFvMesh_H

#include "dynamicFvMesh/dynamicFvMesh.H"
#include "db/dictionary/dictionary.H"
#include "meshes/primitiveShapes/point/pointIOField.H"
#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class multiSliderMotionFvMesh Declaration
\*---------------------------------------------------------------------------*/

class multiSliderMotionFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Dictionary of motion control parameters
        const dictionary dynamicMeshCoeffs_;

        //- Calculated mesh displacement
        vectorField displacement_;

        //- Specified cellZones
        labelList zoneIDs_;

        //- Points to move per cellZone
        labelListList pointIDs_;

        //- moving patches per cellZone
        labelListList movingIDs_;

        //- moving patches velocity or displacment per cellZone
        PtrList<Function1<vector>> motionProfile_;

        //- Whether velocity or displacement profile
        boolList velocityProfile_;

        //- initial displacement per cellZone
        vectorField initialDisplacement_;

        //- buffer region per cellZone
        scalarField bufferSize_;

        //- blend type per cellZone
        wordList blendType_;


public:

    //- Runtime type information
    TypeName("multiSliderMotionFvMesh");


    // Constructors

        //- Construct from IOobject
        multiSliderMotionFvMesh(const IOobject& io);

        //- Disallow default bitwise copy construction
        multiSliderMotionFvMesh(const multiSliderMotionFvMesh&) = delete;


    //- Destructor
    ~multiSliderMotionFvMesh();


    // Member Functions
        autoPtr<indirectPrimitivePatch> makePatch
        (
            const labelList& patchIDs
        );

        void slidePoints
        (
            const vectorField& zoneDisplacements
        );

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiSliderMotionFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
