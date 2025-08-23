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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "rawTopoChangerFvMesh/rawTopoChangerFvMesh.H"
#include "meshes/polyMesh/polyTopoChangeMap/polyTopoChangeMap.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rawTopoChangerFvMesh, 0);
    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        rawTopoChangerFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::rawTopoChangerFvMesh::rawTopoChangerFvMesh(const IOobject& io)
:
    topoChangerFvMesh(io)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rawTopoChangerFvMesh::~rawTopoChangerFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::rawTopoChangerFvMesh::update()
{
    // Do mesh changes (use inflation - put new points in topoChangeMap)
    Info<< "rawTopoChangerFvMesh : Checking for topology changes..."
        << endl;

    // Mesh not moved/changed yet
    moving_ = false;
    topoChanged_ = false;

    // Do any topology changes. Sets topoChanged (through polyTopoChange)
    autoPtr<polyTopoChangeMap> topoChangeMap = topoChanger_.changeMesh(true);

    bool hasChanged = topoChangeMap.valid();

    if (hasChanged)
    {
        Info<< "rawTopoChangerFvMesh : Done topology changes..."
            << endl;

        // Temporary: fix fields on patch faces created out of nothing
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Two situations:
        // - internal faces inflated out of nothing
        // - patch faces created out of previously internal faces

        // Is face mapped in any way?
        PackedBoolList mappedFace(nFaces());

        const label nOldInternal = topoChangeMap().oldPatchStarts()[0];

        const labelList& faceMap = topoChangeMap().faceMap();
        for (label facei = 0; facei < nInternalFaces(); facei++)
        {
            if (faceMap[facei] >= 0)
            {
                mappedFace[facei] = 1;
            }
        }
        for (label facei = nInternalFaces(); facei < nFaces(); facei++)
        {
            if (faceMap[facei] >= 0 && faceMap[facei] >= nOldInternal)
            {
                mappedFace[facei] = 1;
            }
        }

        const List<objectMap>& fromFaces = topoChangeMap().facesFromFacesMap();

        forAll(fromFaces, i)
        {
            mappedFace[fromFaces[i].index()] = 1;
        }

        const List<objectMap>& fromEdges = topoChangeMap().facesFromEdgesMap();

        forAll(fromEdges, i)
        {
            mappedFace[fromEdges[i].index()] = 1;
        }

        const List<objectMap>& fromPts = topoChangeMap().facesFromPointsMap();

        forAll(fromPts, i)
        {
            mappedFace[fromPts[i].index()] = 1;
        }

        // Set unmapped faces to zero
        Info<< "rawTopoChangerFvMesh : zeroing unmapped boundary values."
            << endl;
        zeroUnmappedValues<scalar, fvPatchField, volMesh>(mappedFace);
        zeroUnmappedValues<vector, fvPatchField, volMesh>(mappedFace);
        zeroUnmappedValues<sphericalTensor, fvPatchField, volMesh>(mappedFace);
        zeroUnmappedValues<symmTensor, fvPatchField, volMesh>(mappedFace);
        zeroUnmappedValues<tensor, fvPatchField, volMesh>(mappedFace);

        // Special handling for phi: set unmapped faces to recreated phi
        Info<< "rawTopoChangerFvMesh :"
            << " recreating phi for unmapped boundary values." << endl;
        const volVectorField& U = lookupObject<volVectorField>("U");
        surfaceScalarField& phi = const_cast<surfaceScalarField&>
        (
            lookupObject<surfaceScalarField>("phi")
        );
        setUnmappedValues
        (
            phi,
            mappedFace,
            (linearInterpolate(U) & Sf())()
        );


        if (topoChangeMap().hasMotionPoints())
        {
            pointField newPoints = topoChangeMap().preMotionPoints();

            // Give the meshModifiers opportunity to modify points
            Info<< "rawTopoChangerFvMesh :"
                << " calling modifyMotionPoints." << endl;
            topoChanger_.modifyMotionPoints(newPoints);

            // Actually move points
            Info<< "rawTopoChangerFvMesh :"
                << " calling movePoints." << endl;

            movePoints(newPoints);
        }
    }
    else
    {
        //Pout<< "rawTopoChangerFvMesh :"
        //    << " no topology changes..." << endl;
    }

    return hasChanged;
}


// ************************************************************************* //
