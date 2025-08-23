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

Application

Description

\*---------------------------------------------------------------------------*/

#include "primitives/strings/string/string.H"
#include "db/IOstreams/IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main()
{
    string a("a"), b("b"), c("c"), d("d"), e("e"),
           f("f"), g("g"), h("h"), i("i"), j("j");

    Info<< "1) a = b + c + d + e + f + g;\n";
    a = b + c + d + e + f + g;
    Info<< a << endl;

    {
        Info<< "2) a = e + f + g;\n";
        a = e + f + g;
        Info<< a << endl;
    }

    Info<< "3) a = h + i + j;\n";
    a = h + i + j;
    Info<< a << endl;

    Info<< "end" << endl;

    return 0;
}


// ************************************************************************* //
