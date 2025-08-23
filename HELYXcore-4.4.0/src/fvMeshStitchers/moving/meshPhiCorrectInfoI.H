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

#include "moving/meshPhiCorrectInfo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::meshPhiCorrectInfo::validateFace() const
{
    if (shape_ != shape::face)
    {
        FatalError
            << "Face data requested from a non-face correction info"
            << exit(FatalError);
    }
}


void Foam::meshPhiCorrectInfo::validateCell() const
{
    if (shape_ != shape::cell)
    {
        FatalError
            << "Cell data requested from a non-cell correction info"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshPhiCorrectInfo::meshPhiCorrectInfo()
:
    shape_(shape::invalid),
    delta_(NaN)
{}


Foam::meshPhiCorrectInfo::meshPhiCorrectInfo(const shape& s)
:
    shape_(s),
    delta_(s == shape::invalid ? NaN : 0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::meshPhiCorrectInfo::deltaPhi() const
{
    validateFace();
    return delta_;
}


Foam::scalar& Foam::meshPhiCorrectInfo::deltaPhi()
{
    validateFace();
    return delta_;
}


Foam::scalar Foam::meshPhiCorrectInfo::deltaDVdtError() const
{
    validateCell();
    return delta_;
}


Foam::scalar& Foam::meshPhiCorrectInfo::deltaDVdtError()
{
    validateCell();
    return delta_;
}


template<class TrackingData>
bool Foam::meshPhiCorrectInfo::valid(TrackingData& td) const
{
    return shape_ != shape::invalid;
}


template<class TrackingData>
bool Foam::meshPhiCorrectInfo::sameGeometry
(
    const fvMesh&,
    const meshPhiCorrectInfo& l,
    const scalar tol,
    TrackingData& td
) const
{
    return true;
}


template<class TrackingData>
bool Foam::meshPhiCorrectInfo::updateCell
(
    const fvMesh& mesh,
    const label thisCelli,
    const labelPair& neighbourPatchAndFacei,
    const meshPhiCorrectInfo& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    const label neighbourPatchi = neighbourPatchAndFacei.first();
    const label neighbourFacei = neighbourPatchAndFacei.second();

    const meshPhiPreCorrectInfo& thisPci = td.cellPci_[thisCelli];
    const meshPhiPreCorrectInfo& neighbourPci =
        neighbourPatchi == -1
      ? td.internalFacePci_[neighbourFacei]
      : td.patchFacePci_[neighbourPatchi][neighbourFacei];

    if (!valid(td))
    {
        shape_ = shape::cell;
        deltaDVdtError() = 0;
    }

    if (thisPci.layer() < neighbourPci.layer())
    {
        const label s =
            neighbourPatchi == -1
          ? (mesh.owner()[neighbourFacei] == thisCelli ? +1 : -1)
          : +1;

        deltaDVdtError() -= s*neighbourInfo.deltaPhi();

        return true;
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
bool Foam::meshPhiCorrectInfo::updateFace
(
    const fvMesh& mesh,
    const labelPair& thisPatchAndFacei,
    const label neighbourCelli,
    const meshPhiCorrectInfo& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    const label thisPatchi = thisPatchAndFacei.first();
    const label thisFacei = thisPatchAndFacei.second();

    const meshPhiPreCorrectInfo& thisPci =
        thisPatchi == -1
      ? td.internalFacePci_[thisFacei]
      : td.patchFacePci_[thisPatchi][thisFacei];
    const meshPhiPreCorrectInfo& neighbourPci = td.cellPci_[neighbourCelli];

    if (!valid(td))
    {
        shape_ = shape::face;
        deltaPhi() = 0;
    }

    if (thisPci.layer() < neighbourPci.layer())
    {
        const label s =
            thisPatchi == -1
          ? (mesh.owner()[thisFacei] == neighbourCelli ? +1 : -1)
          : +1;

        deltaPhi() =
            s*thisPci.weight()/neighbourPci.weight()
           *(td.dVdtError_[neighbourCelli] + neighbourInfo.deltaDVdtError());

        return true;
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
bool Foam::meshPhiCorrectInfo::updateFace
(
    const fvMesh& mesh,
    const labelPair& thisPatchAndFacei,
    const meshPhiCorrectInfo& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td))
    {
        shape_ = shape::face;
        deltaPhi() = 0;
    }

    deltaPhi() = - neighbourInfo.deltaPhi();

    return true;
}


template<class TrackingData>
bool Foam::meshPhiCorrectInfo::equal
(
    const meshPhiCorrectInfo& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::meshPhiCorrectInfo::operator==
(
    const Foam::meshPhiCorrectInfo& rhs
) const
{
    return shape_ == rhs.shape_ && delta_ == rhs.delta_;
}


bool Foam::meshPhiCorrectInfo::operator!=
(
    const Foam::meshPhiCorrectInfo& rhs
) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const meshPhiCorrectInfo& l)
{
    return os << static_cast<label>(l.shape_) << token::SPACE << l.delta_;
}


Foam::Istream& Foam::operator>>(Istream& is, meshPhiCorrectInfo& l)
{
    label s;
    is >> s >> l.delta_;
    l.shape_ = static_cast<meshPhiCorrectInfo::shape>(s);
    return is;
}


// ************************************************************************* //
