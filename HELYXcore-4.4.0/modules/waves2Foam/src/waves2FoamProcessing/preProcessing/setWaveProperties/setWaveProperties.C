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
    (c) held by original author

\*---------------------------------------------------------------------------*/

#include "setWaveProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(setWaveProperties, 0);
defineRunTimeSelectionTable(setWaveProperties, setWaveProperties);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void setWaveProperties::lineFormatting( Ostream& os, const word& key)
{
    os << indent << key << token::SPACE;

    for (int i=key.size(); i<Nspaces_-1; i++)
    {
        os << token::SPACE;
    }
}


void setWaveProperties::addITstream
(
    Ostream& os,
    const word& key,
    const ITstream& it
)
{
    lineFormatting(os, key);

    forAll(it, ii)
    {
        os << it[ii];

        if (ii < it.size() - 1)
        {
            os << token::SPACE;
        }
    }

    os << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


void setWaveProperties::writeBeginning( Ostream& os)
{
    wordList names( dict_.name().components('/') );

    fileName ends( names[names.size() -1]);

#if EXTBRANCH==1
    char delim(':');
#elif OFPLUSBRANCH==1
    char delim('.');
#else
    #if OFVERSION<220
        char delim(':');
    #else
        char delim('.');
    #endif
#endif

    wordList subnames( ends.components( delim ) );
    word dictName( subnames[ subnames.size() -1 ] );

    os << dictName << nl << token::BEGIN_BLOCK << incrIndent << nl;
}


void setWaveProperties::writeGiven( Ostream& os )
{
    wordList toc( dict_.toc() );

    forAll(toc, item)
    {
        if (!dict_.isDict(toc[item]))
        {
            ITstream it( dict_.lookup( toc[item] ) );
            addITstream( os, toc[item], it);
        }
    }
}


void setWaveProperties::writeGiven( Ostream& os, word name )
{
    ITstream it( dict_.lookup(name) );

    addITstream( os, name, it );
}


void setWaveProperties::writeDerived( Ostream& os, word name, scalar val)
{
    lineFormatting(os, name);

    os << val << token::END_STATEMENT << nl;
}


void setWaveProperties::writeDerived( Ostream& os, word name, vector val)
{
    lineFormatting(os, name);

    os << val << token::END_STATEMENT << nl;
}


void setWaveProperties::writeDerived( Ostream& os, word name, scalarField val)
{
    lineFormatting(os, name);

    os << "nonuniform List<scalar>" << nl;

    os << indent << val.size() << nl << indent
       << token::BEGIN_LIST << incrIndent << nl;

    forAll(val, vali)
    {
        os << indent << val[vali] << nl;
    }

    os << decrIndent << indent << token::END_LIST
       << token::END_STATEMENT << nl;
}


void setWaveProperties::writeDerived(Ostream& os, word name, vectorField val)
{
    lineFormatting(os, name);

    os << "nonuniform List<vector>" << nl;

    os << indent << val.size() << nl << indent << token::BEGIN_LIST
       << incrIndent << nl;

    forAll(val, vali)
    {
        os << indent << val[vali] << nl;
    }

    os << decrIndent << indent << token::END_LIST
       << token::END_STATEMENT << nl;
}


void setWaveProperties::writeRelaxationZone( Ostream& os )
{
    word rl("relaxationZone");

    if (dict_.found(rl))
    {
        os << nl;
        os.beginBlock(rl);

        dictionary sd(dict_.subDict(rl));

        wordList toc(sd.toc());

        forAll(toc, item)
        {
            ITstream it(sd.lookup(toc[item]));

            addITstream(os, toc[item], it);
        }

        os.endBlock();
    }
}


void setWaveProperties::writeEnding( Ostream& os )
{
    os.endBlock();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


setWaveProperties::setWaveProperties
(
    const Time& rT,
    dictionary& dict,
    bool write
)
:
    write_(write),
    dict_(dict),

    g_
    (
        uniformDimensionedVectorField
        (
            rT.db().lookupObject<uniformDimensionedVectorField>("g")
        ).value()
    )
{
    G_  = Foam::mag(g_);
    PI_ = M_PI;

    Nspaces_ = 20;
}


setWaveProperties::~setWaveProperties()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


autoPtr<setWaveProperties> setWaveProperties::New
(
    const Time& rT,
    dictionary& dict,
    bool write
)
{
    word waveTheoryTypeName = dict.lookup<word>("waveType");

    const auto ctor =
        ctorTableLookup
        (
            "wave property type",
            setWavePropertiesConstructorTable_(),
            waveTheoryTypeName + "Properties"
        );
    return autoPtr<setWaveProperties>(ctor(rT, dict, write));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
