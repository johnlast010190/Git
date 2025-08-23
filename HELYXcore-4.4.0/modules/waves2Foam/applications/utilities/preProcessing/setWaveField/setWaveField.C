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

Application
    setWaveField

Description
    Loop over every cell in the computational domain and set VOF-ratio and
    velocity field accordingly to specified wave theory.

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.

Additional information
    Implementation published and validated in the following journal article:

    @article { jacobsenFuhrmanFredsoe2011,
        Author = {Jacobsen, N G and Fuhrman, D R and Freds\o{}e, J},
        title = {{A Wave Generation Toolbox for the Open-Source CFD Library: OpenFoam\textregistered{}}},
        Journal = {{Int. J. for Numer. Meth. Fluids}},
        Year = {2012},
        Volume = {70},
        Number = {9},
        Pages = {1073-1088},
        DOI = {{10.1002/fld.2726}},
    }

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"

#if EXTBRANCH==1
    #if 310<OFVERSION
        #include "foamTime.H"
    #else
        #include "db/Time/Time.H"
    #endif
#elif OFPLUSBRANCH==1
    #include "db/Time/Time.H"
#else
    #include "db/Time/Time.H"
#endif

//#if EXTBRANCH==1 && OFVERSION>310
//    #include "foamTime.H"
//#else
//    #include "db/Time/Time.H"
//#endif

#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "setWaveField/setWaveField.H"

#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"

#include "include/crossVersionCompatibility.H"

using namespace Foam;

int main(int argc, char *argv[])
{

#include "include/addTimeOptions.H"
#include "include/setRootCase.H"

#include "include/createTime.H"
#include "include/createMesh.H"

#include "cfdTools/general/include/readGravitationalAcceleration.H"

#include "include/readWaveProperties.H"

    Info<< "\nReading field alpha\n" << endl;
    volScalarField alpha
    (
        IOobject
        (
            Foam::waves2Foam::aName(),
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field p\n" << endl;
    volScalarField pd
    (
        IOobject
        (
            Foam::waves2Foam::pName(),
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Setting the wave field ...\n" << endl;

    setWaveField swf(mesh, U, alpha, pd);

    swf.correct();

    alpha.write();

    U.write();

    pd.write();

    Info<< nl << "End" << nl << endl;

    return 0;
}
