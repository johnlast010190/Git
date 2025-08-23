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
    (c) 2020 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "reaction/reaction/reaction.H"
#include "db/dictionary/dictionary.H"


// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

Foam::label Foam::reaction::nUnNamedReactions(0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::reaction::getNewReactionID()
{
    return nUnNamedReactions++;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reaction::reaction
(
    const speciesTable& species,
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs
)
:
    name_("un-named-reaction-" + Foam::name(getNewReactionID())),
    species_(species),
    pyrolisisGases_(),
    lhs_(lhs),
    rhs_(rhs),
    glhs_(),
    grhs_(),
    gaseousSpcecies_(false)
{}


Foam::reaction::reaction
(
    const reaction& r,
    const speciesTable& species
)
:
    name_(r.name() + "Copy"),
    species_(species),
    pyrolisisGases_(r.gasSpecies()),
    lhs_(r.lhs_),
    rhs_(r.rhs_),
    glhs_(r.glhs_),
    grhs_(r.grhs_),
    gaseousSpcecies_(r.gaseousSpcecies_)
{}


Foam::reaction::reaction
(
    const speciesTable& species,
    const dictionary& dict
)
:
    name_(dict.dictName()),
    species_(species)
{
    specieCoeffs::setLRhs
    (
        IStringStream(dict.lookup("reaction"))(),
        species_,
        lhs_,
        rhs_
    );

    // Solid reaction
    gaseousSpcecies_ = dict.parent().parent().found("gaseousSpecies");
    if (gaseousSpcecies_)
    {
        pyrolisisGases_ = dict.parent().parent().lookup("gaseousSpecies");
        specieCoeffs::setLRhs
        (
            IStringStream(dict.lookup("reaction"))(),
            pyrolisisGases_,
            glhs_,
            grhs_
        );
    }
}


void Foam::reaction::solidReactionStrLeft(OStringStream& reaction) const
{
    for (label i = 0; i < glhs().size(); ++i)
    {
        if (i > 0)
        {
            reaction << " + ";
        }
        if (mag(glhs()[i].stoichCoeff - 1) > SMALL)
        {
            reaction << glhs()[i].stoichCoeff;
        }
        reaction << gasSpecies()[glhs()[i].index];
        if
        (
            mag(glhs()[i].exponent.operator scalar() - glhs()[i].stoichCoeff)
          > SMALL
        )
        {
            reaction << "^" << glhs()[i].exponent;
        }
    }
}


void Foam::reaction::solidReactionStrRight(OStringStream& reaction) const
{
    for (label i = 0; i < grhs().size(); ++i)
    {
        if (i > 0)
        {
            reaction << " + ";
        }
        if (mag(grhs()[i].stoichCoeff - 1) > SMALL)
        {
            reaction << grhs()[i].stoichCoeff;
        }
        reaction << gasSpecies()[grhs()[i].index];
        if
        (
            mag(grhs()[i].exponent.operator scalar() - grhs()[i].stoichCoeff)
          > SMALL
        )
        {
            reaction << "^" << grhs()[i].exponent;
        }
    }
}


Foam::string Foam::reaction::solidReactionStr(OStringStream& reaction) const
{
    specieCoeffs::reactionStr(reaction, this->species(), this->lhs());
    if (glhs().size() > 0)
    {
        reaction << " + ";
        solidReactionStrLeft(reaction);
    }
    reaction << " = ";
    specieCoeffs::reactionStr(reaction, this->species(), this->rhs());
    if (grhs().size() > 0)
    {
        reaction << " + ";
        solidReactionStrRight(reaction);
    }
    return reaction.str();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::reaction::write(Ostream& os) const
{
    OStringStream reaction;

    if (gasSpecies().size())
    {
        os.writeEntry("reaction", solidReactionStr(reaction));
    }
    else
    {
        os.writeEntry
        (
            "reaction",
            specieCoeffs::reactionStr(reaction, species_, lhs_, rhs_)
        );
    }
}


// ************************************************************************* //