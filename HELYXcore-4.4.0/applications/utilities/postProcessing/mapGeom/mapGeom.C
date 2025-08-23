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
    Find the best map of flow configuration from the database to the
    current case based on machine-learning method.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "sample/helyxSample/helyxSample.H"
#include "triSurface/triSurface.H"
#include "containers/Array3d/Array3d.H"
#include "helyxMap.H"
#include "sample/trainGeoModel/trainGeoModel.H"

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addNote("Map geometric object using surface mesh file in constant/triSurface");

    #include "include/setRootCase.H"

    Foam::Info<< "Create time\n" << Foam::endl;

	Foam::Time runTime(Foam::Time::controlDictName, args);


	IOdictionary dict
	(
		IOobject
		(
            "sampleGeomDict",
			runTime.system(),
			"",
			runTime,
			IOobject::MUST_READ,
			IOobject::NO_WRITE,
			false
		)
	);


	Info<< "Time now = " << runTime.timeName() << endl;

    word home="~/";
    fileName databaseDir(home / "helyxAI"/"database"/"car3d");

    fileName dbsDir = dict.lookupOrDefault<fileName> ("databaseDir","default");
    if (dbsDir !="default")
    {
        databaseDir=dbsDir;
    }

    databaseDir= databaseDir.expand();
    Info<<"databaseDir:"<<databaseDir<<endl;

    fileName caseListFile(databaseDir / "caselist.dbs");
    fileName geomlistFile(databaseDir / "objectMask.dbs");
    bool exists = isFile(geomlistFile.c_str());

	if (!exists)
	{
		FatalErrorInFunction<<"Error mapping flow field, database file: "
		<<geomlistFile<<" not found or cannot open."
		<< exit(FatalError);
	}

    trainGeoModel train
    (
		&runTime,
		dict
	);


    Info<<"Read database, this will take a while...."<<endl;
	train.readDatabase(geomlistFile);

    vector ndiv0(127, 63, 63);

    vector ndiv=dict.lookupOrDefault<vector> ("ndivisions",ndiv0);

    label nx=label(ndiv[0]);
    label ny=label(ndiv[1]);
    label nz=label(ndiv[2]);

    Info<<"Grid divisions: "<<nx<<" "<<ny<<" "<<nz<<endl;

    train.buildModel(nx, ny, nz);

    train.getPrediction();
	Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	  << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	  << nl << endl;

    Info<< "End\n" << endl;

	return 0;
}

// ************************************************************************* //
