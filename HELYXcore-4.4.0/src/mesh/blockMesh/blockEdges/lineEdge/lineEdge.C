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
    (c) 2019-2021 OpenCFD Ltd

\*---------------------------------------------------------------------------*/

#include "blockEdges/lineEdge/lineEdge.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blockEdges
{
    defineTypeNameAndDebug(lineEdge, 0);
    addToRunTimeSelectionTable(blockEdge, lineEdge, Istream);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockEdges::lineEdge::lineEdge
(
    const pointField& points,
    const edge& fromTo
)
:
    blockEdge(points, fromTo)
{}


Foam::blockEdges::lineEdge::lineEdge
(
    const pointField& points,
    const label from,
    const label to
)
:
    blockEdge(points, from, to)
{}


Foam::blockEdges::lineEdge::lineEdge
(
    const dictionary& dict,
    const label index,
    const searchableSurfaces&,
    const pointField& points,
    Istream& is
)
:
    blockEdge(dict, index, points, is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::blockEdges::lineEdge::position(const scalar lambda) const
{
    return blockEdge::linearPosition(lambda);
}


Foam::scalar Foam::blockEdges::lineEdge::length() const
{
    return Foam::mag(lastPoint() - firstPoint());
}


// ************************************************************************* //
