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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::label Foam::pairPotentialList::pairPotentialIndex
(
    const label a,
    const label b
) const
{
    label index;

    if (a < b)
    {
        index = a*(2*nIds_ - a - 1)/2 + b;
    }

    else
    {
        index = b*(2*nIds_ - b - 1)/2 + a;
    }

    if (index > size() - 1)
    {
        FatalErrorInFunction
            << "Attempting to access a pairPotential with too high an index."
            << nl << "a = " << a << ", b = " << b << ", index = " << index
            << nl << "max index = " << size() - 1
            << nl << abort(FatalError);
    }

    return index;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::pairPotentialList::rCutMax() const
{
    return rCutMax_;
}


inline Foam::scalar Foam::pairPotentialList::rCutMaxSqr() const
{
    return rCutMaxSqr_;
}


inline const Foam::pairPotential& Foam::pairPotentialList::electrostatic() const
{
    return electrostaticPotential_;
}


// ************************************************************************* //
