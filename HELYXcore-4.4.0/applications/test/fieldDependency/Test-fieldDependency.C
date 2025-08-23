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
    Test field dependencies.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "fields/volFields/volFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Creating field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0)
    );

    Info<< "Creating field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0)
    );


    Info<< "p.eventNo:" << p.eventNo() << endl;
    Info<< "p.uptodate:" << p.upToDate(T)<< endl;

    // Change T and mark as uptodate.
    Info<< "Changing T" << endl;
    T = 0.0;
    T.setUpToDate();
    Info<< "T.eventNo:" << T.eventNo() << endl;

    // Check p dependency:
    Info<< "p.uptodate:" << p.upToDate(T)<< endl;

    // Change p and mark as uptodate.
    Info<< "Changing p." << endl;
    p.setUpToDate();
    Info<< "p.uptodate:" << p.upToDate(T)<< endl;
    Info<< "p.eventNo:" << p.eventNo() << endl;


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
