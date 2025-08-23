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
    (c) 2019-2022 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "reaction/specieCoeffs/specieCoeffs.H"
#include "containers/Lists/DynamicList/DynamicList.H"
#include "db/IOstreams/StringStreams/OStringStream.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::specieCoeffs::specieCoeffs
(
    const speciesTable& species,
    Istream& is
)
{
    token t(is);
    if (t.isNumber())
    {
        stoichCoeff = t.number();
        is >> t;
    }
    else
    {
        stoichCoeff = 1;
    }

    exponent = stoichCoeff;

    if (t.isWord())
    {
        word specieName = t.wordToken();

        size_t i = specieName.find('^');

        if (i != word::npos)
        {
            string exponentStr = specieName
            (
                i + 1,
                specieName.size() - i - 1
            );
            exponent = scalar(atof(exponentStr.c_str()));
            specieName = specieName(0, i);
        }

        if (species.found(specieName))
        {
            index = species[specieName];
        }
        else
        {
            // In order to support the convoluted way in which
            // the solidChemistry currently handles gaseous species
            // an error cannot be generated here.
            // We will re-introduce this check after solidChemistry is
            // re-written in a rational manner.
            // FatalIOErrorInFunction(is)
            //     << "Specie " << specieName
            //     << " not found in table " << species
            //     << exit(FatalIOError);

            index = -1;
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Expected a word but found " << t.info()
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::specieCoeffs::setLRhs
(
    Istream& is,
    const speciesTable& species,
    List<specieCoeffs>& lhs,
    List<specieCoeffs>& rhs
)
{
    DynamicList<specieCoeffs> dlrhs;

    while (is.good())
    {
        dlrhs.append(specieCoeffs(species, is));

        if (dlrhs.last().index != -1)
        {
            token t(is);
            if (t.isPunctuation())
            {
                if (t == token::ADD)
                {
                }
                else if (t == token::ASSIGN)
                {
                    lhs = dlrhs.shrink();
                    dlrhs.clear();
                }
                else
                {
                    rhs = dlrhs.shrink();
                    is.putBack(t);
                    return;
                }
            }
            else
            {
                rhs = dlrhs.shrink();
                is.putBack(t);
                return;
            }
        }
        else
        {
            dlrhs.remove();
            if (is.good())
            {
                token t(is);
                if (t.isPunctuation())
                {
                    if (t == token::ADD)
                    {
                    }
                    else if (t == token::ASSIGN)
                    {
                        lhs = dlrhs.shrink();
                        dlrhs.clear();
                    }
                    else
                    {
                        rhs = dlrhs.shrink();
                        is.putBack(t);
                        return;
                    }
                }
            }
            else
            {
                if (!dlrhs.empty())
                {
                    rhs = dlrhs.shrink();
                }
                return;
            }
        }
    }

    FatalIOErrorInFunction(is)
        << "Cannot continue reading reaction data from stream"
        << exit(FatalIOError);
}


void Foam::specieCoeffs::reactionStr
(
    OStringStream& reaction,
    const speciesTable& species,
    const List<specieCoeffs>& scs
)
{
    for (label i = 0; i < scs.size(); ++i)
    {
        if (i > 0)
        {
            reaction << " + ";
        }
        if (mag(scs[i].stoichCoeff - 1) > SMALL)
        {
            reaction << scs[i].stoichCoeff;
        }
        reaction << species[scs[i].index];
        if (mag(scs[i].exponent.operator scalar() - scs[i].stoichCoeff) > SMALL)
        {
            reaction << "^" << scs[i].exponent;
        }
    }
}


Foam::string Foam::specieCoeffs::reactionStr
(
    OStringStream& reaction,
    const speciesTable& species,
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs
)
{
    reactionStr(reaction, species, lhs);
    reaction << " = ";
    reactionStr(reaction, species, rhs);
    return reaction.str();
}


// ************************************************************************* //