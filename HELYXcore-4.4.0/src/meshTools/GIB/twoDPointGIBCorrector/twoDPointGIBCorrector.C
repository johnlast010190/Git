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
    (c) 2017 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "GIB/twoDPointGIBCorrector/twoDPointGIBCorrector.H"
#include "meshes/polyMesh/polyMesh.H"
#include "meshes/polyMesh/polyPatches/constraint/wedge/wedgePolyPatch.H"
#include "meshTools/meshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twoDPointGIBCorrector::twoDPointGIBCorrector(const polyMesh& mesh)
:
    twoDPointCorrector(mesh)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void twoDPointGIBCorrector::correctPoints
(
    pointField& p,
    const pointField& baseP
) const
{
    if (!required()) return;
    // Algorithm:
    // Loop through all edges. Calculate the average point position A for
    // the front and the back. Correct the position of point P (in two planes)
    // such that vectors AP and planeNormal are parallel

    // Get reference to edges
    const edgeList&  meshEdges = mesh().edges();

    const labelList& neIndices = normalEdgeIndices();
    const vector& pn = planeNormal();
    vector pnn = cmptMultiply(pn, pn);

    forAll(neIndices, edgeI)
    {
        const label& lStart = meshEdges[neIndices[edgeI]].start();
        const label& lEnd = meshEdges[neIndices[edgeI]].end();
        point& pStart = p[lStart];
        point& pEnd = p[lEnd];

        // calculate average point position
        point A = 0.5*(pStart + pEnd);
        meshTools::constrainToMeshCentre(mesh_, A);

        if (isWedge())
        {
            snapToWedge(pn, A, pStart);
            snapToWedge(pn, A, pEnd);
        }
        else
        {
            pStart -= cmptMultiply(pStart, pnn);
            pEnd = pStart;

            pStart += cmptMultiply(baseP[lStart], pnn);
            pEnd += cmptMultiply(baseP[lEnd], pnn);
        }
    }
}


void twoDPointGIBCorrector::correctPlanes(pointField& p) const
{
    if (!required()) return;

    // Get reference to edges
    const edgeList&  meshEdges = mesh().edges();

    const labelList& neIndices = normalEdgeIndices();
    const vector& pn = planeNormal();
    vector pnn = cmptMultiply(pn, pn);

    forAll(neIndices, edgeI)
    {
        const label& lStart = meshEdges[neIndices[edgeI]].start();
        const label& lEnd = meshEdges[neIndices[edgeI]].end();
        point& pStart = p[lStart];
        point& pEnd = p[lEnd];

        vector pStartCor = cmptMultiply(pStart, pnn);
        pStart -= pStartCor;
        pEnd = pStart;
    }
}


void twoDPointGIBCorrector::correctPointsPlus(pointField& p) const
{
    twoDPointCorrector::correctPoints(p);
    const polyBoundaryMesh& patches = mesh().boundaryMesh();
    forAll(patches, patchI)
    {
        if (isA<wedgePolyPatch>(patches[patchI]))
        {
            vector pn = refCast<const wedgePolyPatch>(patches[patchI]).n();
            const labelList& wpp = patches[patchI].meshPoints();
            forAll(wpp, pI)
            {
                label gpI = wpp[pI];
                vector disp =  p[gpI] - mesh().points()[gpI];
                vector dispV= (disp&pn)*pn;
                disp -= dispV;
                p[gpI] = mesh().points()[gpI] + disp;
            }
        }
    }
}


void twoDPointGIBCorrector::correctPointsPlus
(
    pointField& p,
    const pointField& baseP
) const
{
    correctPoints(p, baseP);
    const polyBoundaryMesh& patches = mesh().boundaryMesh();
    forAll(patches, patchI)
    {
        if (isA<wedgePolyPatch>(patches[patchI]))
        {
            vector pn = refCast<const wedgePolyPatch>(patches[patchI]).n();
            const labelList& wpp = patches[patchI].meshPoints();
            forAll(wpp, pI)
            {
                label gpI = wpp[pI];
                vector disp =  p[gpI] - mesh().points()[gpI];
                vector dispV= (disp&pn)*pn;
                disp -= dispV;
                p[gpI] = mesh().points()[gpI] + disp;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
