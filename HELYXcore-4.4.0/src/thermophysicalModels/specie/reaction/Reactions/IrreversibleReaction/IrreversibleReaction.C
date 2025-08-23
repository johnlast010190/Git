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
    (c) 2011-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "reaction/Reactions/IrreversibleReaction/IrreversibleReaction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType, class ReactionRate>
Foam::IrreversibleReaction<ThermoType, ReactionRate>::IrreversibleReaction
(
    const Reaction<ThermoType>& reaction,
    const ReactionRate& k
)
:
    Reaction<ThermoType>(reaction),
    k_(k)
{}


template<class ThermoType, class ReactionRate>
Foam::IrreversibleReaction<ThermoType, ReactionRate>::IrreversibleReaction
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
Foam::IrreversibleReaction<ThermoType, ReactionRate>::IrreversibleReaction
(
    const IrreversibleReaction<ThermoType, ReactionRate>& irr,
    const speciesTable& species
)
:
    Reaction<ThermoType>(irr, species),
    k_(irr.k_)
{}


template<class ThermoType, class ReactionRate>
Foam::IrreversibleReaction<ThermoType, ReactionRate>::IrreversibleReaction
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
void Foam::IrreversibleReaction<ThermoType, ReactionRate>::preEvaluate() const
{
    k_.preEvaluate();
}


template<class ThermoType, class ReactionRate>
void Foam::IrreversibleReaction<ThermoType, ReactionRate>::postEvaluate() const
{
    k_.postEvaluate();
}


template<class ThermoType, class ReactionRate>
Foam::scalar Foam::IrreversibleReaction<ThermoType, ReactionRate>::kf
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
Foam::scalar Foam::IrreversibleReaction<ThermoType, ReactionRate>::kr
(
    const scalar kfwd,
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return 0;
}


template<class ThermoType, class ReactionRate>
Foam::scalar Foam::IrreversibleReaction<ThermoType, ReactionRate>::kr
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return 0;
}


template<class ThermoType, class ReactionRate>
Foam::scalar Foam::IrreversibleReaction<ThermoType, ReactionRate>::dkfdT
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
Foam::scalar Foam::IrreversibleReaction<ThermoType, ReactionRate>::dkrdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    const scalar dkfdT,
    const scalar kr
) const
{
    return 0;
}


template<class ThermoType, class ReactionRate>
bool Foam::IrreversibleReaction<ThermoType, ReactionRate>::hasDkdc() const
{
    return k_.hasDdc();
}


template<class ThermoType, class ReactionRate>
void Foam::IrreversibleReaction<ThermoType, ReactionRate>::dkfdc
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
void Foam::IrreversibleReaction<ThermoType, ReactionRate>::dkrdc
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
    dkrdc = 0;
}



template<class ThermoType, class ReactionRate>
void Foam::IrreversibleReaction<ThermoType, ReactionRate>::write
(
    Ostream& os
) const
{
    Reaction<ThermoType>::write(os);
    k_.write(os);
}


// ************************************************************************* //
