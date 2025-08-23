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
    (c) 2011-2022 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "reaction/Reactions/ReversibleReaction/ReversibleReaction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType, class ReactionRate>
Foam::ReversibleReaction<ThermoType, ReactionRate>::ReversibleReaction
(
    const Reaction<ThermoType>& reaction,
    const ReactionRate& k
)
:
    Reaction<ThermoType>(reaction),
    k_(k)
{}


template<class ThermoType, class ReactionRate>
Foam::ReversibleReaction<ThermoType, ReactionRate>::ReversibleReaction
(
    const speciesTable& species,
    const PtrList<ThermoType>& speciesThermo,
    const dictionary& dict
)
:
    Reaction<ThermoType>(species, speciesThermo, dict),
    k_(species, dict)
{}


template<class ThermoType, class ReactionRate>
Foam::ReversibleReaction<ThermoType, ReactionRate>::ReversibleReaction
(
    const ReversibleReaction<ThermoType, ReactionRate>& rr,
    const speciesTable& species
)
:
    Reaction<ThermoType>(rr, species),
    k_(rr.k_)
{}


template<class ThermoType, class ReactionRate>
Foam::ReversibleReaction<ThermoType, ReactionRate>::ReversibleReaction
(
    const speciesTable& species,
    const PtrList<ThermoType>& speciesThermo,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    Reaction<ThermoType>(species, speciesThermo, dict),
    k_(species, obr, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType, class ReactionRate>
void Foam::ReversibleReaction<ThermoType, ReactionRate>::preEvaluate() const
{
    k_.preEvaluate();
}


template<class ThermoType, class ReactionRate>
void Foam::ReversibleReaction<ThermoType, ReactionRate>::postEvaluate() const
{
    k_.postEvaluate();
}


template<class ThermoType, class ReactionRate>
Foam::scalar Foam::ReversibleReaction<ThermoType, ReactionRate>::kf
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return k_(p, T, c, li);
}


template<class ThermoType, class ReactionRate>
Foam::scalar Foam::ReversibleReaction<ThermoType, ReactionRate>::kr
(
    const scalar kfwd,
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return kfwd/max(this->Kc(p, T), ROOTSMALL);
}


template<class ThermoType, class ReactionRate>
Foam::scalar Foam::ReversibleReaction<ThermoType, ReactionRate>::kr
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return kr(kf(p, T, c, li), p, T, c, li);
}

template<class ThermoType, class ReactionRate>
Foam::scalar Foam::ReversibleReaction<ThermoType, ReactionRate>::dkfdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return k_.ddT(p, T, c, li);
}


template<class ThermoType, class ReactionRate>
Foam::scalar Foam::ReversibleReaction<ThermoType, ReactionRate>::dkrdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    const scalar dkfdT,
    const scalar kr
) const
{
    const scalar Kc = max(this->Kc(p, T), ROOTSMALL);

    return dkfdT/Kc - (Kc > ROOTSMALL ? kr*this->dKcdTbyKc(p, T) : 0);
}


template<class ThermoType, class ReactionRate>
bool Foam::ReversibleReaction<ThermoType, ReactionRate>::hasDkdc() const
{
    return k_.hasDdc();
}


template<class ThermoType, class ReactionRate>
void Foam::ReversibleReaction<ThermoType, ReactionRate>::dkfdc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dkfdc
) const
{
    k_.ddc(p, T, c, li, dkfdc);
}


template<class ThermoType, class ReactionRate>
void Foam::ReversibleReaction<ThermoType, ReactionRate>::dkrdc
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
    const scalar Kc = max(this->Kc(p, T), ROOTSMALL);
    dkrdc = dkfdc/Kc;
}


template<class ThermoType, class ReactionRate>
void Foam::ReversibleReaction<ThermoType, ReactionRate>::write
(
    Ostream& os
) const
{
    Reaction<ThermoType>::write(os);
    k_.write(os);
}


// ************************************************************************* //
