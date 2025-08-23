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
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::potential::nIds() const
{
    return idList_.size();
}


inline const Foam::List<Foam::word>& Foam::potential::idList() const
{
    return idList_;
}


inline const Foam::List<Foam::word>& Foam::potential::siteIdList() const
{
    return siteIdList_;
}


inline Foam::scalar Foam::potential::potentialEnergyLimit() const
{
    return potentialEnergyLimit_;
}


inline Foam::label Foam::potential::nPairPotentials() const
{
    return pairPotentials_.size();
}


inline const Foam::labelList& Foam::potential::removalOrder() const
{
    return removalOrder_;
}


inline const Foam::pairPotentialList& Foam::potential::pairPotentials() const
{
    return pairPotentials_;
}


inline const Foam::tetherPotentialList&
Foam::potential::tetherPotentials() const
{
    return tetherPotentials_;
}


inline const Foam::vector& Foam::potential::gravity() const
{
    return gravity_;
}


// ************************************************************************* //
