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
    (c) 2022 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "moving/meshPhiPreCorrectInfo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshPhiPreCorrectInfo::meshPhiPreCorrectInfo()
:
    layer_(-labelMax),
    weight_(NaN)
{}


Foam::meshPhiPreCorrectInfo::meshPhiPreCorrectInfo
(
    const label layer,
    const scalar weight
)
:
    layer_(layer),
    weight_(weight)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::meshPhiPreCorrectInfo::layer() const
{
    return layer_;
}


Foam::scalar Foam::meshPhiPreCorrectInfo::weight() const
{
    return weight_;
}


template<class TrackingData>
bool Foam::meshPhiPreCorrectInfo::valid(TrackingData& td) const
{
    return layer_ != -labelMax;
}


template<class TrackingData>
bool Foam::meshPhiPreCorrectInfo::sameGeometry
(
    const fvMesh&,
    const meshPhiPreCorrectInfo& l,
    const scalar tol,
    TrackingData& td
) const
{
    return layer_ == l.layer_;
}


template<class TrackingData>
bool Foam::meshPhiPreCorrectInfo::updateCell
(
    const fvMesh& mesh,
    const label thisCelli,
    const labelPair& neighbourPatchAndFacei,
    const meshPhiPreCorrectInfo& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td))
    {
        layer_ = neighbourInfo.layer_ + 1;
        weight_ = 0;
    }

    if (layer_ > neighbourInfo.layer_)
    {
        weight_ += neighbourInfo.weight_;

        return true;
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
bool Foam::meshPhiPreCorrectInfo::updateFace
(
    const fvMesh& mesh,
    const labelPair& thisPatchAndFacei,
    const label neighbourCelli,
    const meshPhiPreCorrectInfo& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td))
    {
        layer_ = neighbourInfo.layer_ + 1;
    }

    if (layer_ > neighbourInfo.layer_)
    {
        weight_ = neighbourInfo.weight_;

        return true;
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
bool Foam::meshPhiPreCorrectInfo::updateFace
(
    const fvMesh& mesh,
    const labelPair& thisPatchAndFacei,
    const meshPhiPreCorrectInfo& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (layer_ < neighbourInfo.layer_)
    {
        layer_ = neighbourInfo.layer_;
        weight_ = neighbourInfo.weight_;
    }

    return true;
}


template<class TrackingData>
bool Foam::meshPhiPreCorrectInfo::equal
(
    const meshPhiPreCorrectInfo& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::meshPhiPreCorrectInfo::operator==
(
    const Foam::meshPhiPreCorrectInfo& rhs
) const
{
    return layer_ == rhs.layer_ && weight_ == rhs.weight_;
}


bool Foam::meshPhiPreCorrectInfo::operator!=
(
    const Foam::meshPhiPreCorrectInfo& rhs
) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const meshPhiPreCorrectInfo& l)
{
    return os << l.layer_ << token::SPACE << l.weight_;
}


Foam::Istream& Foam::operator>>(Istream& is, meshPhiPreCorrectInfo& l)
{
    return is >> l.layer_ >> l.weight_;
}


// ************************************************************************* //
