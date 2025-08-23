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

#include "hessianMeshOptimization/layerCellHandling/layerPointStack/layerPointStack.H"
#include "db/IOstreams/IOstreams/Istream.H"
#include "db/IOstreams/IOstreams/Ostream.H"
#include "db/IOstreams/token/token.H"
#include "primitives/contiguous/contiguous.H"

namespace Foam
{

    layerPointStack::layerPointStack(Istream& is)
    {
        operator>>(is, *this);
    }

    Istream& operator>>(Istream& is, layerPointStack& fI)
    {
        bool bs;
        List<vector> sP;
        labelList pointS;
        label lNu;

        is >> bs >> sP >> pointS >> lNu;

        fI = layerPointStack(bs, sP, pointS, lNu);
        return is;
    }

    Ostream& operator<<(Ostream& os, const layerPointStack& fI)
    {
        os     <<fI.baseStack() <<" "
            <<fI.stackPoints()<< " "
            <<fI.pointStackLabels()<<" "
            <<fI.layerNu()<<" ";

        return os;
    }
}
