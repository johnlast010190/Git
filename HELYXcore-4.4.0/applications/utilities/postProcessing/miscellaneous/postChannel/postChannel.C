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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Application
    postChannel

Group
    grpPostProcessingUtilities

Description
    Post-processes data from channel flow calculations.

    For each time: calculate: txx, txy, tyy, txy,
    eps, prod, vorticity, enstrophy and helicity. Assuming that the mesh
    is periodic in the x and z directions, collapse Umeanx, Umeany, txx,
    txy and tyy to a line and print them as standard output.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "channelIndex.H"
#include "graphField/makeGraph.H"
#include "include/OSspecific.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    timeSelector::addOptions();

    #include "include/setRootCase.H"
    Info<< "Create time\n" << endl;
    Time runTime(Time::controlDictName, args);

    // Get times list
    instantList timeDirs = timeSelector::select0(runTime, args);

    Info<< "Create mesh for time = " << runTime.timeName() << nl << endl;
    fvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
           IOobject::MUST_READ
        )
    );

    const word& gFormat = runTime.graphFormat();

    // Setup channel indexing for averaging over channel down to a line
    IOdictionary channelDict
    (
        IOobject
        (
            "postChannelDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    channelIndex channelIndexing(mesh, channelDict);


    // For each time step read all fields
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Collapsing fields for time " << runTime.timeName() << endl;

        IOobject UMeanHeader
        (
            "UMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        if (!UMeanHeader.typeHeaderOk<volVectorField>(true))
        {
            Info<< "    No UMean field" << endl;
            continue;
        }

        volVectorField UMean(UMeanHeader, mesh);

        volSymmTensorField UPrime2Mean
        (
            IOobject
            (
                "UPrime2Mean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );
        volScalarField Rxx(UPrime2Mean.component(symmTensor::XX));
        volScalarField Ryy(UPrime2Mean.component(symmTensor::YY));
        volScalarField Rzz(UPrime2Mean.component(symmTensor::ZZ));
        volScalarField Rxy(UPrime2Mean.component(symmTensor::XY));

        volScalarField pPrime2Mean
        (
            IOobject
            (
                "pPrime2Mean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );

        // Average fields over channel down to a line
        fileName path
        (
            UMean.rootPath()/UMean.caseName()/"graphs"/UMean.instance()
        );
        mkDir(path);

        scalarField UMeanXvalues
        (
            channelIndexing.collapse(UMean.component(vector::X)())
        );

        scalarField UMeanYvalues
        (
            channelIndexing.collapse(UMean.component(vector::Y)())
        );

        scalarField UMeanZvalues
        (
            channelIndexing.collapse(UMean.component(vector::Z)())
        );

        scalarField RxxValues(channelIndexing.collapse(Rxx));
        scalarField RyyValues(channelIndexing.collapse(Ryy));
        scalarField RzzValues(channelIndexing.collapse(Rzz));
        scalarField RxyValues(channelIndexing.collapse(Rxy, true));

        scalarField pPrime2MeanValues(channelIndexing.collapse(pPrime2Mean));

        scalarField urmsValues(sqrt(mag(RxxValues)));
        scalarField vrmsValues(sqrt(mag(RyyValues)));
        scalarField wrmsValues(sqrt(mag(RzzValues)));

        scalarField kValues
        (
            0.5*(sqr(urmsValues) + sqr(vrmsValues) + sqr(wrmsValues))
        );


        const scalarField& y = channelIndexing.y();

        makeGraph(y, UMeanXvalues, "Uf", path, gFormat);
        makeGraph(y, urmsValues, "u", path, gFormat);
        makeGraph(y, vrmsValues, "v", path, gFormat);
        makeGraph(y, wrmsValues, "w", path, gFormat);
        makeGraph(y, RxyValues, "uv", path, gFormat);
        makeGraph(y, kValues, "k", path, gFormat);

        makeGraph(y, pPrime2MeanValues, "pPrime2Mean", path, gFormat);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
