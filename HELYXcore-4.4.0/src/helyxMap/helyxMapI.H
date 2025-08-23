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
    (c) 2019-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::helyxMap::searchNearest
(
    std::vector<label>& nb,
    std::vector<scalar>& dist,
    scalar wdisti,
    scalar x,
    scalar y,
    scalar z
)
{
    if (!wDistMap_)
    {
        return searchNearest(nb, dist, x, y, z);
    }

    unsigned int ic = label(wdisti/deltaY_);
    std::vector<scalar> qvect(3);
    qvect[0] = x;
    qvect[1] = y;
    qvect[2] = z;

    if (ic > kdTrees_.size() - 1)
    {
        ic = kdTrees_.size() - 1;
    }

    if (!kdTrees_[ic].active)
    {
        return internalMap_.search(qvect, nb, dist);
    }

    if (wdisti > sourceAlphaMax_ || ic > kdTrees_.size() - 1)
    {
        return internalMap_.search(qvect, nb, dist);
    }
    else
    {
        if (kdTrees_[ic].active)
        {
            kdTrees_[ic].search(qvect, nb, dist);
        }
        else
        {
            internalMap_.search(qvect, nb, dist);
        }
    }

    return;
}


inline void Foam::helyxMap::searchNearest
(
    std::vector<label>& nb,
    std::vector<scalar>& dist,
    scalar x,
    scalar y,
    scalar z
)
{
    std::vector<scalar> qvect(3);
    qvect[0] = x;
    qvect[1] = y;
    qvect[2] = z;

    internalMap_.search(qvect, nb, dist);
}


// ************************************************************************* //
