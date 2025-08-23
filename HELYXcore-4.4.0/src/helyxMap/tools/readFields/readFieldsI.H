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

inline Foam::label Foam::readFields::cellIndex(scalar alphai, label kdTreeId)
{
    if (alphai > alphaMax_)
    {
        return cellMap_[kdTreeId];
    }
    else
    {
        label iAlpha = int(alphai/deltaY_);

        if (kdTreeId < label(alphaMap_[iAlpha].size()))
        {
            return alphaMap_(iAlpha, kdTreeId);
        }
        else
        {
            return -1;
        }
    }
}


inline void Foam::readFields::searchNearest
(
    std::vector<label>& nb,
    std::vector<scalar>& dist,
    scalar alphai,
    scalar x,
    scalar y,
    scalar z
)
{
    if (type_ != "source") return;

    std::vector<scalar> qvect(3);
    qvect[0] = x;
    qvect[1] = y;
    qvect[2] = z;

    if (alphai > alphaMax_)
    {
        internalMap_.search(qvect, nb, dist);
    }
    else
    {
        label ic = label(alphai/(deltaY_ + VSMALL));

        kdTrees_[ic].search(qvect, nb, dist);
    }
}


inline void Foam::readFields::searchNearestBFace
(
    std::vector<label>& nb,
    std::vector<scalar>& dist,
    scalar x,
    scalar y,
    scalar z
)
{
    if (type_ != "source") return;

    std::vector<scalar> qvect(3);
    qvect[0] = x;
    qvect[1] = y;
    qvect[2] = z;

    boundaryMap_.search(qvect, nb, dist);
}


inline bool Foam::readFields::inDomain(const point& pt)
{
    return bBox().contains(pt);
}


// ************************************************************************* //
