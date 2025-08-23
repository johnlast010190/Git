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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/primitiveMesh/PatchTools/PatchTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class FaceList, class PointField>
Foam::labelList Foam::PatchTools::edgeOwner
(
    const PrimitivePatch<FaceList, PointField>& p
)
{
    typedef typename PrimitivePatch<FaceList, PointField>::FaceType FaceType;

    const edgeList& edges = p.edges();
    const labelListList& edgeFaces = p.edgeFaces();
    const List<FaceType>& localFaces = p.localFaces();

    // create the owner list
    labelList edgeOwner(edges.size(), -1);

    forAll(edges, edgeI)
    {
        const labelList& nbrFaces = edgeFaces[edgeI];

        if (nbrFaces.size() == 1)
        {
            edgeOwner[edgeI] = nbrFaces[0];
        }
        else
        {
            // Find the first face whose vertices are aligned with the edge.
            // with multiply connected edges, this is the best we can do
            forAll(nbrFaces, i)
            {
                const FaceType& f = localFaces[nbrFaces[i]];

                if (f.edgeDirection(edges[edgeI]) > 0)
                {
                    edgeOwner[edgeI] = nbrFaces[i];
                    break;
                }
            }

            if (edgeOwner[edgeI] == -1)
            {
                FatalErrorInFunction
                    << "Edge " << edgeI << " vertices:" << edges[edgeI]
                    << " is used by faces " << nbrFaces
                    << " vertices:"
                    << UIndirectList<FaceType>(localFaces, nbrFaces)()
                    << " none of which use the edge vertices in the same order"
                    << nl << "I give up" << abort(FatalError);
            }
        }
    }

    return edgeOwner;
}


// ************************************************************************* //
