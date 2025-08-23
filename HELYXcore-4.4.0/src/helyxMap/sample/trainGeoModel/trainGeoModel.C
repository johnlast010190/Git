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
    (c) 2009-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "sample/trainGeoModel/trainGeoModel.H"
#include "sample/npGrid/npGrid.H"
#include "triSurface/triSurface.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::trainGeoModel::trainGeoModel()
:
    helyxSample()
{}


Foam::trainGeoModel::trainGeoModel(const Time* runTime, const dictionary& dict)
:
    helyxSample(runTime, nullptr)
{
    word gFile = dict.lookupOrDefault<word>("geofile", "car.stl");

    word cdir = cwd();
    caseDir_ = cdir;
    fileName bfile(cdir/"constant"/"triSurface"/gFile);

    triSurface surfaces(bfile);

    const vectorField& cf = surfaces.Cf();

    Info<< "Number of surface elements: " << cf.size() << endl;

    point pmin = min(cf);
    point pmax = max(cf);
    point dp = pmax - pmin;

    Info<< "Bounding box for surface mesh: " << pmin << " " << pmax << endl;
    Info<< "Span of geometry: " << dp << endl;

    point p0(-0.5, -0.5, -0.1);
    point p1(4, 2.2, 2);

    boundBox b0(p0, p1);
    boundBox bbox = dict.lookupOrDefault<boundBox>("boundBox", b0);

    Info<< "Boundbox for surface grid: " << bbox << endl;

    vector ndiv0(127, 63, 63);
    vector ndiv = dict.lookupOrDefault<vector>("ndivisions", ndiv0);

    label nx = label(ndiv[0]);
    label ny = label(ndiv[1]);
    label nz = label(ndiv[2]);

    Info<< "Grid divisions: " << nx << " " << ny << " " << nz << endl;

    Array3d<label> mask;
    mask.setSize(nx, ny, nz);
    mask = 0;

    pmin = bbox.min();
    pmax = bbox.max();
    dp = pmax - pmin;

    scalar dx = dp[0]/float(nx - 1);
    scalar dy = dp[1]/float(ny - 1);

    label nzpt = max(1, nz - 1);
    scalar dz = dp[2]/float(nzpt);

    // Update mask
    forAll(cf, ip)
    {
        point pt = cf[ip];
        if (!bbox.contains(pt))
        {
            continue;
        }

        label i = label((pt[0] - pmin[0])/dx);
        label j = label((pt[1] - pmin[1])/dy);
        label k = label((pt[2] - pmin[2])/dz);

        i = min(i, nx - 1);
        j = min(j, ny - 1);
        k = min(k, nz - 1);

        mask(i, j, k) = 1;
    }

    label cnt = 0;
    for (label i = 0; i < nx; i++)
    {
        for (label j = 0; j < ny; j++)
        {
            for (label k = 0; k < nz; k++)
            {
                cnt += mask(i, j, k);
            }
        }
    }

    gridField grid(bbox, nx - 1, ny - 1, nzpt);

    grid.setFieldSize();
    grid.inSolids_ = mask;

    sampleGrids_["geomObject"] = grid;

    Info<< "Number of surface grids: " << cnt << endl;
    Info<< "Total grids: " << nx*ny*nz << endl;
}


Foam::trainGeoModel::trainGeoModel(const Time* runTime, const fvMesh* mesh)
:
    helyxSample(runTime, mesh)
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::trainGeoModel::~trainGeoModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::trainGeoModel::readDatabase1(const word& dbsfile)
{
	word input;
	npGrid::readTextFile(input, dbsfile);

	wordList vtmp;
	stringTokenize(vtmp, input, "\n");

	vCaseItems_.resize(0);
	for (label i = 0; i < vtmp.size() - 1; i += 2)
	{
		const word& s1 = vtmp[i];
		const word& s2 = vtmp[i + 1];

		if (s1.length() < 3 || s2.size() < 3) break;

		dataItem ditem;
		ditem.setDataInfo(s1);
		ditem.setData(s2);

		vCaseItems_.push_back(ditem);
	}

	Info<< "Number of database items read: " << vCaseItems_.size() << endl;
}


void Foam::trainGeoModel::readDatabase(const word& dbsfile)
{
	std::ifstream is(dbsfile.c_str());
	if (!is)
	{
		FatalErrorInFunction
			<< "Error reading database, database file: "
			<< dbsfile << " not found or cannot be opened."
			<< exit(FatalError);
	}

	word line;
	word cnts;

	vCaseItems_.resize(0);
	while (!is.eof())
	{
		std::getline(is, line, '\n');
		if (line.length() < 3) break;

		dataItem ditem;
		ditem.setDataInfo(line);

		std::getline(is, cnts, '\n');
		if (cnts.length() < 3) break;

		ditem.setData(cnts);
		vCaseItems_.push_back(ditem);
	}

	Info<< "Number of database items read: " << vCaseItems_.size() << endl;
}


void Foam::trainGeoModel::setOptions
(
	const DynamicList<word>& argNames,
	const DynamicList<word>& argValues
)
{
	if (argNames.size() == 0)
    {
        return;
    }

	for (label i = 0; i < argNames.size(); i++)
	{
		word key = argNames[i];
		word val = argValues[i];

		if (key == "-databaseDir")
		{
			helyxSample::databaseDir_ = val;

			continue;
		}

		if (key == "-task") continue;
	}
}


void Foam::trainGeoModel::getCaseLabels()
{
	caseNames_.resize(0);

    std::map<word, dataItem>::iterator it;
	for (it = caseItems_.begin(); it != caseItems_.end(); it++)
	{
		caseNames_.push_back(it->first);
	}

	std::sort(caseNames_.begin(), caseNames_.end());
	for (unsigned i = 0; i < caseNames_.size(); i++)
	{
		caseLabels_[caseNames_[i]] = i;
	}

	for (it = caseItems_.begin(); it != caseItems_.end(); it++)
	{
		word cname = it->first;
		it->second.identity(caseLabels_[cname]);
	}
}


void Foam::trainGeoModel::readCaseList()
{
	word caseFile = helyxSample::databaseDir_ + "/caselist.dbs";

	std::ifstream is(caseFile.c_str());
	if (!is)
	{
		FatalErrorInFunction
			<< "Error reading training case list, \n"
            << caseFile << " cannot be opened for reading."
			<< exit(FatalError);
	}

	word str;
	while (!is.eof())
	{
		std::getline(is, str, '\n');
		if (str.length() < 6) continue;

		dataItem ditem(str);
		caseItems_[ditem.caseName()] = ditem;
	}

	is.close();

	Info<< "Number of cases read: " << caseItems_.size() << endl;
}


void Foam::trainGeoModel::readVector
(
	std::vector<scalar>& vect,
    std::ifstream& is
)
{
	word str;

	for (unsigned i = 0; i < vect.size(); i++)
	{
		is  >> vect[i];
	}

    // Eat the end of line sign
	is  >> str;
}


void Foam::trainGeoModel::clearGeoData()
{
	for (unsigned i = 0; i < vCaseItems_.size(); i++)
	{
		vCaseItems_[i].clearGeoData();
	}
}


void Foam::trainGeoModel::buildModel(label m, label n, label k)
{
	model_ = new KNN();

	farray2d trainData;

	label cnt = 0;
	for (unsigned i = 0; i < vCaseItems_.size(); i++)
	{
		if
		(
			vCaseItems_[i].m() != m
		 || vCaseItems_[i].n() != n
		 || vCaseItems_[i].k() != k
		)
		{
			continue;
		}

		trainData.push_back(vCaseItems_[i].geoData_);

		idmap_[cnt] = i;
		cnt++;
	}

	if (trainData.empty())
	{
		Info<< "Warning: no data available to build KNN search tree." << endl;
	}
	else
	{
		model_().build(trainData);
	}

	Info<< "Number of training data: " << trainData.size() << endl;
}


void Foam::trainGeoModel::buildModel()
{
	model_ = new KNN();

	label ncases = vCaseItems_.size();
	farray2d trainData;

	std::map<word, dataItem>::iterator it;
	for (label i = 0; i < ncases; i++)
	{
		word cname = caseNames_[i];
		if (vCaseItems_[i].geoData_.size() > 10)
		{
			trainData.push_back(vCaseItems_[i].geoData_);
		}
	}

	if (trainData.empty())
	{
		Info<< "Warning: no data available to build KNN search tree." << endl;
	}
	else
	{
		model_().build(trainData);
	}

	clearGeoData();
}


void Foam::trainGeoModel::readGeometryFile()
{
	word geofile = helyxSample::databaseDir_ + "/geom.dbs";

	std::ifstream is(geofile.c_str());
	if (!is)
	{
		FatalErrorInFunction
			<< "Error reading geometry representation file, \n"
            << geofile << " cannot be opened for reading."
			<< exit(FatalError);
	}

	caseNames_.resize(0);
	vCaseItems_.resize(0);

	word caseName;
	word sDimension;
	word sData;

	Info<< "Reading geometric representation file, this will take some time..."
        << endl;

	while (!is.eof())
	{
		std::getline(is, caseName, '\n');
		if (caseName.size() < 3) break;

		std::getline(is, sDimension, '\n');
		if (sDimension.size() < 3) break;

		std::getline(is, sData, '\n');
		if (sData.size() < 3) break;

		wordList vtmp;
		helyxMap::stringTokenize(vtmp, sDimension, " ");

		if (vtmp.size() != 4)
		{
			FatalErrorInFunction
                << "Error reading geometry representation file, "
                << "dimension data line " << sDimension << " is incorrect"
                << exit(FatalError);
		}

		label npoints = atoi(vtmp[0].c_str());

		label m = atoi(vtmp[1].c_str()) + 1;
		label n = atoi(vtmp[2].c_str()) + 1;
		label k = atoi(vtmp[3].c_str());

		if (k > 1)
		{
			k += 1;
		}

		std::vector<scalar> vdata(npoints);

		wordList svdata;

		split(svdata, sData, ' ');

		label nsample = svdata.size();
		if (nsample != npoints)
		{
			FatalErrorInFunction
                << "Error reading geometry representation file, "
                << "incorrect data item for " << caseName
                << exit(FatalError);
		}

		for (label i = 0; i < npoints; i++)
		{
			vdata[i] = atof(svdata[i].c_str());
		}

		dataItem ditem;
		ditem.setData(m, n, k, vdata);

		vCaseItems_.push_back(ditem);
		caseNames_.push_back(caseName);
	}

	is.close();

	Info<< "Number of cases read: " << vCaseItems_.size() << endl;
}


void Foam::trainGeoModel::saveModel()
{
	fileName geomdl = helyxSample::databaseDir_ + "/geoModel.mdl";

	Info<< "Save geometry recognition model into database..." << endl;

	model_().saveModel(geomdl);
}


Foam::scalar Foam::trainGeoModel::getDist(label i, label j)
{
	const std::vector<scalar> v1 = vCaseItems_[i].geoData_;
	const std::vector<scalar> v2 = vCaseItems_[j].geoData_;

	scalar sum = 0;

	for (unsigned t = 0; t < v1.size(); t++)
	{
		sum += fabs(v1[t] - v2[t]);
	}

	return sum;
}


void Foam::trainGeoModel::getGeomSample()
{
	std::map<word, gridField>::iterator it;

	for (it = sampleGrids_.begin(); it != sampleGrids_.end(); it++)
	{
		helyxSample::getGeomSample(it->second);
	}
}


void Foam::trainGeoModel::getPrediction()
{
	std::map<word, gridField>::iterator it;

	for (it = sampleGrids_.begin(); it != sampleGrids_.end(); it++)
	{
		gridField& grid = it->second;
		word cnts = "mask";

		for (label i = 0; i < grid.nx_ + 1; ++i)
		{
			for (label j = 0; j < grid.ny_ + 1; ++j)
			{
				for (label k = 0; k < grid.nzpt(); ++k)
				{
					label inSolid = grid.inSolids_(i, j, k);

					cnts += ("," + std::to_string(inSolid));
				}
			}
		}

		wordList vtmp0;
		stringTokenize(vtmp0, cnts, ",");

		std::vector<scalar> vv0(vtmp0.size() - 1);
		for (label s = 1; s < vtmp0.size(); ++s)
		{
			vv0[s - 1] = atof(vtmp0[s].c_str());
		}

		std::vector<int> nb(2);
		std::vector<scalar> vdists(2);
		std::vector<scalar> vquery = vv0;

		model_().search(vquery, nb, vdists);

		label id = nb[0];
		label id1 = nb[1];

		if (idmap_.size() == vCaseItems_.size())
		{
			id = idmap_[id];
			id1 = idmap_[id1];

			Info<< "Using idmap: " << id << " " << id1 << endl;
		}

		vCaseItems_[id].show();
		Info<< "Distance: " << vdists[0] << nl << endl;

		vCaseItems_[id1].show();
		Info<< "Distance: " << vdists[1] << nl << endl;
	}
}


void Foam::trainGeoModel::getSolidLabels
(
	gridField& grid,
	std::vector<label>& inSolid
)
{
	Info<< "Total number of grid nodes: " << grid.size() << endl;

	inSolid.resize(grid.size());

	label cnt = 0;
	for (label i = 0; i < grid.nx_ + 1; ++i)
	{
		for (label j = 0; j < grid.ny_ + 1; ++j)
		{
			for (label k = 0; k < grid.nzpt(); ++k)
			{
				label idx = grid.ijk(i, j, k);
				inSolid[idx] = grid.inSolids_(i, j, k);

				cnt += inSolid[idx];
			}
		}
	}

	Info<< "Insolids: " << cnt << endl;
}


void Foam::trainGeoModel::loadModel()
{
	Info<< "Loading pre-trained flow-topology-recognition model..." << endl;

	word geomdl = helyxSample::databaseDir_ + "/geoModel.mdl";

	model_ = new KNN(geomdl);
}


// ************************************************************************* //
