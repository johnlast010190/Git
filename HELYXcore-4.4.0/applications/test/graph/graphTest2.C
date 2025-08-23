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
    (c) 2011-2017 OpenFOAM Foundation

Application
    graphTest

Description
    Test program for making graphs

\*---------------------------------------------------------------------------*/

#include "graph.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main()
{
    scalarField eta(200);
    scalarField C1mR(eta.size());

    forAll(eta, i)
    {
        eta[i] = scalar(i)/10.0;
    }

    scalar C1 = 1.42;
    scalar eta0 = 4.38;
    scalar beta = 0.012;

    C1mR = C1 - ((eta*(1.0 - eta/eta0))/(1.0 + beta*pow(eta, 3.0)));

    graph("C&1! - R", "C&1! - R", "@h!", eta,  C1mR).write
    (
        "C1mR",
        "xmgr"
    );

    Info<< "end" << endl;
}


// ************************************************************************* //
