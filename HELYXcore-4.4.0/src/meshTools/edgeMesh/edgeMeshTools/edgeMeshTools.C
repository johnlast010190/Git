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
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "edgeMesh/edgeMeshTools/edgeMeshTools.H"

#include "edgeMesh/extendedEdgeMesh/extendedFeatureEdgeMesh/extendedFeatureEdgeMesh.H"
#include "db/IOstreams/Fstreams/OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::edgeMeshTools::writeStats
(
    Ostream& os,
    const extendedFeatureEdgeMesh& emesh
)
{
    os  << "Feature set:" << nl
        << "    points : " << emesh.points().size() << nl
        << "    of which" << nl
        << "        convex             : "
        << emesh.concaveStart() << nl
        << "        concave            : "
        << (emesh.mixedStart()-emesh.concaveStart()) << nl
        << "        mixed              : "
        << (emesh.nonFeatureStart()-emesh.mixedStart()) << nl
        << "        non-feature        : "
        << (emesh.points().size()-emesh.nonFeatureStart()) << nl
        << "    edges  : " << emesh.edges().size() << nl
        << "    of which" << nl
        << "        external edges     : "
        << emesh.internalStart() << nl
        << "        internal edges     : "
        << (emesh.flatStart()- emesh.internalStart()) << nl
        << "        flat edges         : "
        << (emesh.openStart()- emesh.flatStart()) << nl
        << "        open edges         : "
        << (emesh.multipleStart()- emesh.openStart()) << nl
        << "        multiply connected : "
        << (emesh.edges().size()- emesh.multipleStart()) << endl;
}


// ************************************************************************* //
