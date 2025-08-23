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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2022-2023 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "reaction/Reactions/ReactionList/ReactionList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::ReactionList<ThermoType>::ReactionList
(
    const speciesTable& species,
    const PtrList<ThermoType>& speciesThermo,
    const dictionary& dict
)
{
    // Set general temperature limits from the dictionary
    Reaction<ThermoType>::TlowDefault =
        dict.lookupOrDefault<scalar>("Tlow", 0);

    Reaction<ThermoType>::ThighDefault =
        dict.lookupOrDefault<scalar>("Thigh", GREAT);

    if (dict.found("reactions"))
    {
        const dictionary& reactions = dict.subDict("reactions");

        this->setSize(reactions.size());
        label i = 0;

        forAllConstIter(dictionary, reactions, iter)
        {
            this->set
            (
                i++,
                Reaction<ThermoType>::New
                (
                    species,
                    speciesThermo,
                    reactions.subDict(iter().keyword())
                ).ptr()
            );
        }
    }
}


template<class ThermoType>
Foam::ReactionList<ThermoType>::ReactionList
(
    const speciesTable& species,
    const PtrList<ThermoType>& speciesThermo,
    const objectRegistry& obr,
    const dictionary& dict
)
{
    // Set general temperature limits from the dictionary
    Reaction<ThermoType>::TlowDefault =
        dict.lookupOrDefault<scalar>("Tlow", 0);

    Reaction<ThermoType>::ThighDefault =
        dict.lookupOrDefault<scalar>("Thigh", GREAT);

    const dictionary& reactions(dict.subDict("reactions"));

    this->setSize(reactions.size());
    label i = 0;

    forAllConstIter(dictionary, reactions, iter)
    {
        this->set
        (
            i++,
            Reaction<ThermoType>::New
            (
                species,
                speciesThermo,
                obr,
                reactions.subDict(iter().keyword())
            ).ptr()
        );
    }
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::ReactionList<ThermoType>::write(Ostream& os) const
{
    if ((*this).size())
    {
        os.beginBlock("reactions");
        forAll(*this, i)
        {
            const Reaction<ThermoType>& r = this->operator[](i);
            os.beginBlock(r.name());
            os.writeEntry("type", r.type());
            r.write(os);
            os.endBlock();
        }
        os.endBlock();
    }

    os.writeEntry("Tlow", Reaction<ThermoType>::TlowDefault);
    os.writeEntry("Thigh", Reaction<ThermoType>::ThighDefault);

    os << nl;
}


// ************************************************************************* //
