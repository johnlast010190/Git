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
    (c) 2018-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "ReactionProxy.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::ReactionProxy<ThermoType>::ReactionProxy
(
    const speciesTable& species,
    const PtrList<ThermoType>& speciesThermo,
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs
)
:
    Reaction<ThermoType>(species, speciesThermo, lhs, rhs)
{}


template<class ThermoType>
Foam::autoPtr<Foam::Reaction<ThermoType>>
Foam::ReactionProxy<ThermoType>::clone() const
{
    NotImplemented;
}


template<class ThermoType>
Foam::autoPtr<Foam::Reaction<ThermoType>>
Foam::ReactionProxy<ThermoType>::clone(const speciesTable& species) const
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::ReactionProxy<ThermoType>::preEvaluate() const
{}


template<class ThermoType>
void Foam::ReactionProxy<ThermoType>::postEvaluate() const
{}


template<class ThermoType>
Foam::scalar Foam::ReactionProxy<ThermoType>::kf
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    NotImplemented;
}


template<class ThermoType>
Foam::scalar Foam::ReactionProxy<ThermoType>::kr
(
    const scalar kfwd,
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    NotImplemented;
}


template<class ThermoType>
Foam::scalar Foam::ReactionProxy<ThermoType>::kr
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    NotImplemented;
}


template<class ThermoType>
Foam::scalar Foam::ReactionProxy<ThermoType>::dkfdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    NotImplemented;
}


template<class ThermoType>
Foam::scalar Foam::ReactionProxy<ThermoType>::dkrdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    const scalar dkfdT,
    const scalar kr
) const
{
    NotImplemented;
}


template<class ThermoType>
bool Foam::ReactionProxy<ThermoType>::hasDkdc() const
{
    NotImplemented;
}


template<class ThermoType>
void Foam::ReactionProxy<ThermoType>::dkfdc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dkfdc
) const
{
    NotImplemented;
}


template<class ThermoType>
void Foam::ReactionProxy<ThermoType>::dkrdc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    const scalarField& dkfdc,
    const scalar kr,
    scalarField& dkrdc
) const
{
    NotImplemented;
}


// ************************************************************************* //