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

Description

\*---------------------------------------------------------------------------*/

#include "primitives/strings/lists/stringListOps.H"
#include "db/IOstreams/StringStreams/IStringStream.H"
#include "db/IOstreams/IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    stringList strLst
    (
        IStringStream
        (
            "("
            "\"hello\""
            "\"heello\""
            "\"heeello\""
            "\"bye\""
            "\"bbye\""
            "\"bbbye\""
            "\"okey\""
            "\"okkey\""
            "\"okkkey\""
            ")"
        )()
    );

    wordReList reLst(IStringStream("( okey \"[hy]e+.*\" )")());

    Info<< "stringList " << strLst << nl;

    labelList matches = findStrings(".*ee.*", strLst);

    Info<< "matches found for regexp .*ee.* :" << nl << matches << nl;
    forAll(matches, i)
    {
        Info<< " -> " << strLst[matches[i]] << nl;
    }
    Info<< endl;

    matches = findStrings(reLst, strLst);

    Info<< "matches found for " << reLst << nl << matches << nl;
    forAll(matches, i)
    {
        Info<< " -> " << strLst[matches[i]] << nl;
    }
    Info<< endl;

    stringList subLst = subsetStrings(".*ee.*", strLst);
    Info<< "subset stringList: " << subLst << nl;

    subLst = subsetStrings(reLst, strLst);
    Info<< "subset stringList: " << subLst << nl;

    inplaceSubsetStrings(reLst, strLst);
    Info<< "subsetted stringList: " << strLst << nl;

    inplaceSubsetStrings(".*l.*", strLst);
    Info<< "subsetted stringList: " << strLst << nl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
