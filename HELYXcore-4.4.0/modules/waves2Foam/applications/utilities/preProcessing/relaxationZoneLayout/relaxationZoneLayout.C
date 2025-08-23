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
    relaxationZoneLayout

Description
    Loop over every cell in the computational domain and set the relaxation zone
    coordinate sigma and show which cells belong to which relaxatio zone.

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


#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "relaxationZone/relaxationShape/relaxationShape.H"
#include "relaxationZone/relaxationWeight/relaxationWeight.H"

#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"

using namespace Foam;

int main(int argc, char *argv[])
{

#include "include/addTimeOptions.H"
#include "include/setRootCase.H"

#include "include/createTime.H"

#include "include/createMesh.H"

#include "cfdTools/general/include/readGravitationalAcceleration.H"

#include "include/readWaveProperties.H"

    Info<< "Creating field relaxationZoneLayout\n" << endl;

    volScalarField rzl
    (
        IOobject
        (
            "relaxationZoneLayout",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, -1),
        "zeroGradient"
    );

    Info<< "Creating field relaxationZoneLayoutSigma\n" << endl;

    volScalarField srzl
    (
        IOobject
        (
            "relaxationZoneSigmaValue",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0),
        "zeroGradient"
    );

    volScalarField wrzl
    (
        IOobject
        (
            "relaxationZoneWeightOnComputed",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 1),
        "zeroGradient"
    );

    // Read the relaxation names from waveProperties
    wordList relaxNames(waveProperties.lookup("relaxationNames"));

    // Create the relaxation shape function and add the data to the above
    // generated fields
    forAll(relaxNames, relaxi)
    {
        // Get relaxation zone cells and coorsponding sigma coordinates
        autoPtr<relaxationShapes::relaxationShape> relaxShape
            = relaxationShapes::relaxationShape::New(relaxNames[relaxi], mesh);

        const labelList& cells(relaxShape->cells());

        const scalarField& sigma(relaxShape->sigma());

        forAll(cells, celli)
        {
            rzl[cells[celli]] += (relaxi + 1);
            srzl[cells[celli]] = sigma[celli];
        }

        Info<< relaxNames[relaxi] << " has "
             << cells.size() << " cells\n" << endl;

        // Compute relaxation zone weights
        autoPtr<relaxationWeights::relaxationWeight> relaxWeight
            = relaxationWeights::relaxationWeight::
              New(relaxNames[relaxi], mesh);

        scalarField weights(sigma.size(), 1.0);

        relaxWeight->weights(cells, sigma, weights);

        forAll(weights, celli)
        {
            wrzl[cells[celli]] = weights[celli];
        }
    }

    Info<< "Write the fields to Time = " << runTime.timeName() << endl;
    rzl.write();
    srzl.write();
    wrzl.write();

    Info<< nl << "End" << endl;

    return 0;
}
