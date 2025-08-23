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
    (c) 2019-2021 Engys Ltd.

Description
    Sample a block around a solid object with uniform grid to generate
    the mask of the object.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "helyxMap.H"
#include "sample/helyxSample/helyxSample.H"
#include "sample/helyxSamplePar/helyxSamplePar.H"

int main(int argc, char *argv[])
{

	#include "include/setRootCase.H"

    Foam::Info<< "Create time\n" << Foam::endl;

	Foam::Time runTime(Foam::Time::controlDictName, args);


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

    IOdictionary dict
	(
		IOobject
		(
			"helyxMapDict",
			runTime.system(),
			"",
			runTime,
			IOobject::MUST_READ,
			IOobject::NO_WRITE,
			false
		)
	);


    helyxSample *gridSample;

    if (!Pstream::parRun())
    {
        gridSample =
            new helyxSample
            (
                &runTime,
                &mesh
            );
    }
    else
    {
        gridSample =
            new helyxSamplePar
            (
                &runTime,
                &mesh
            );
    }

    gridSample->setInput(dict, args);
    gridSample->nProc_=Pstream::nProcs();
    gridSample->caseDir_=cwd();

    if (Pstream::parRun())
    {
        gridSample->myId_     = Pstream::myProcNo();
        gridSample->masterId_ = Pstream::masterNo();
        gridSample->parInit(Pstream::nProcs());
    }


    instantList sourceTimes = runTime.times();
    label sourceTimeIndex = runTime.timeIndex();

    if (gridSample->mapTimeName_ == "latestTime")
    {
        sourceTimeIndex = sourceTimes.size() - 1;
    }
    else
    {
        IStringStream is(gridSample->mapTimeName_);
        sourceTimeIndex = Time::findClosestTimeIndex
        (
            sourceTimes,
            readScalar(is)
        );
    }

    runTime.setTime(sourceTimes[sourceTimeIndex], sourceTimeIndex);
    gridSample->mapTimeName_ = runTime.timeName();

    Info<< "\nSampling time " << runTime.timeName() << endl;

    gridSample->source_ =
        new readFields
        (
            &mesh,
            &runTime,
            gridSample
        );

    //gridSample->source().createFields(gridSample->mapTimeName_);
    gridSample->source().setInputs();

    if (!Pstream::parRun())
    {
        gridSample->buildSearchTrees();
        gridSample->getSampleData("geometry");
    }
    else
    {
        // const boundBox &bbox=gridSample->sampleBox();
        // gridSample->source().storeFields(bbox);

        gridSample->wait();

       // gridSample->combineFields();
        gridSample->getBodySurface();

        gridSample->wait();
        Info<<"construct knn...."<<endl;
        if (gridSample->master())
        {
            gridSample->constructKnn();
            gridSample->getSampleData("geometry");
        }

    }

    if (Pstream::parRun())
    {
        gridSample->wait();
    }


    if (!Pstream::parRun() || gridSample->master())
    {

        gridSample->saveObjectMasks(gridSample->caseType_);
    }

    Info<<"Geometry mask saved to file: "<<gridSample->databasePath()<<endl;
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	  << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	  << nl << endl;

    Info<< "End\n" << endl;

    delete gridSample;

	return 0;
}

// ************************************************************************* //
