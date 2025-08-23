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
    (c) 2019 Engys Ltd.
    (c) 2014-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "primitives/Pair/PairKey/PairKey.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PairKey<Type>::hash::hash()
{}


template<class Type>
Foam::PairKey<Type>::PairKey(const word& orderedKeyword)
:
    ordered_(false),
    orderedKeyword_(orderedKeyword)
{}


template<class Type>
Foam::PairKey<Type>::PairKey
(
    const word& name1,
    const word& name2,
    const bool ordered,
    const word& orderedKeyword
)
:
    Pair<word>(name1, name2),
    ordered_(ordered),
    orderedKeyword_(orderedKeyword)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::PairKey<Type>::~PairKey()
{}


// * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::PairKey<Type>::read(Istream& is, const bool exitOnError)
{
    token nextToken(is);
    if
    (
        !is.good()
     || !nextToken.isPunctuation()
     || nextToken.pToken() != token::BEGIN_LIST
    )
    {
        if (exitOnError)
        {
            FatalIOErrorInFunction(is)
                << "Expected '" << token::BEGIN_LIST
                << "' when reading pair type."
                << exit(FatalIOError);
        }
        else
        {
            return false;
        }
    }

    word f, s;
    bool ordered;

    if (exitOnError)
    {
        is >> f;
    }
    else
    {
        is >> nextToken;
        if (is.good() && nextToken.isWord())
        {
            f = nextToken.wordToken();
        }
        else
        {
            return false;
        }
    }

    is >> nextToken;
    if (!is.good())
    {
        if (exitOnError)
        {
            FatalIOErrorInFunction(is)
                << "Expected pair of words when reading pair type."
                << exit(FatalIOError);
        }
        else
        {
            return false;
        }
    }
    if (nextToken.isWord() && nextToken.wordToken() == orderedKeyword_)
    {
        ordered = true;
    }
    else if (nextToken.isWord() && nextToken.wordToken() == "and")
    {
        ordered = false;
    }
    else
    {
        ordered = false;
        is.putBack(nextToken);
    }

    if (exitOnError)
    {
        is >> s;
    }
    else
    {
        is >> nextToken;
        if (is.good() && nextToken.isWord())
        {
            s = nextToken.wordToken();
        }
        else
        {
            return false;
        }
    }

    is >> nextToken;
    if
    (
        !is.good()
     || !(nextToken.isPunctuation() && nextToken.pToken() == token::END_LIST)
    )
    {
        if (exitOnError)
        {
            FatalIOErrorInFunction(is)
                << "Pair type is not recognised. "
                << "Use (name1 " << orderedKeyword_
                << " name2) for an ordered pair, and "
                << "(name1 name2) or (name1 and name2) for an unordered pair."
                << exit(FatalIOError);
        }
        else
        {
            return false;
        }
    }

    this->first() = f;
    this->second() = s;
    ordered_ = ordered;

    return true;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
Foam::label Foam::PairKey<Type>::hash::operator()
(
    const PairKey<Type>& key
) const
{
    if (key.ordered_)
    {
        return
            typename Type::hash()
            (
                key.first(),
                typename Type::hash()(key.second())
            );
    }
    else
    {
        return
            typename Type::hash()(key.first())
          + typename Type::hash()(key.second());
    }
}


// * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * * //

template<class Type>
bool Foam::operator==
(
    const PairKey<Type>& a,
    const PairKey<Type>& b
)
{
    const label c = Pair<Type>::compare(a,b);

    return
        (a.ordered_ == b.ordered_)
     && (
            (a.ordered_ && (c == 1))
         || (!a.ordered_ && (c != 0))
        );
}


template<class Type>
bool Foam::operator!=
(
    const PairKey<Type>& a,
    const PairKey<Type>& b
)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * * //

template<class Type>
Foam::Istream& Foam::operator>>(Istream& is, PairKey<Type>& key)
{
    key.read(is, true);
    return is;
}


// * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const PairKey<Type>& key)
{
    os  << token::BEGIN_LIST
        << key.first()
        << token::SPACE;
    if (key.ordered_)
    {
        os  << key.orderedKeyword_
            << token::SPACE;
    }
    os  << key.second()
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
