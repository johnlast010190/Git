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
    (c) 2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "blocks/namedBlock/namedBlock.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blocks
{
    defineTypeNameAndDebug(namedBlock, 0);
    addToRunTimeSelectionTable(block, namedBlock, Istream);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blocks::namedBlock::namedBlock
(
    const dictionary& dict,
    const label index,
    const pointField& vertices,
    const blockEdgeList& edges,
    const blockFaceList& faces,
    Istream& is
)
:
    word(is),
    block(dict, index, vertices, edges, faces, is)
{
    dictionary& d = const_cast<dictionary&>(dict);
    dictionary* varDictPtr = d.subDictPtr("namedBlocks");
    if (varDictPtr)
    {
        const_cast<dictionary&>(*varDictPtr).add(*this, index);
    }
    else
    {
        dictionary varDict;
        varDict.add(*this, index);
        d.add("namedBlocks", varDict);
    }
}


// ************************************************************************* //
