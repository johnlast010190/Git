/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.0.1
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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#pragma once

// OpenFOAM includes
#include "primitives/Vector/Vector.H"
#include "types/blockMeshBlockType.H"

namespace Foam::functionObjects::runTimeVis
{

struct BlockMeshBlock
{
    BlockMeshBlockType type;
    List<label> points;
    Vector<label> elements;
    word style;
    List<scalar> expansionRatios;

    BlockMeshBlock() = default;
    explicit BlockMeshBlock(Istream& is)
    {
        read(is);
    }

    void read(Istream& is)
    {
        is >> type;
        is >> points;
        is >> elements;
        is >> style;
        is >> expansionRatios;

        // Check state of Istream
        is.check(FUNCTION_NAME);
    }
};

Istream& operator>>(Istream& is, BlockMeshBlock& m);

} // End namespace

// ************************************************************************* //
