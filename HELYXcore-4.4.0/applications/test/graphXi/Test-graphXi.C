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
    graphTest

Description
    Test program for making graphs

\*---------------------------------------------------------------------------*/

#include "graph/graph.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "global/constants/mathematical/mathematicalConstants.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main()
{
    scalarField x(100);
    scalarField r(x.size());

    forAll(x, i)
    {
        x[i] = -3 + 0.06*i;
    }

    scalarField b(0.5*(1.0 + erf(x)));
    scalarField c(1.0 - b);
    scalarField gradb((1/Foam::sqrt(constant::mathematical::pi))*exp(-sqr(x)));
    scalarField lapb(-2*x*gradb);

    r = lapb*b*c/(gradb*gradb);


    graph("r", "x", "r", x, r).write("r", "xmgr");

    Info<< "end" << endl;

    return 0;
}


// ************************************************************************* //
