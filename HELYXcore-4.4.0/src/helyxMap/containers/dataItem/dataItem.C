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

#include "containers/dataItem/dataItem.H"
#include "helyxMap.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dataItem::dataItem()
:
	caseName_("block1"),
	rhoRef_(1.205),
	caseType_("2dcar"),
	identity_(0)
{
	Uref_.x() = 30;
	Uref_.y() = 0;
	Uref_.z() = 0;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::dataItem::getDist(const std::vector<scalar>& vdata) const
{
	scalar sum = 0;
	for (unsigned i = 0; i < vdata.size(); i++)
	{
		sum += fabs(vdata[i] - geoData_[i]);
	}

	return sum;
}


Foam::label Foam::dataItem::numSolids() const
{
	scalar sum = 0;
	for (unsigned i = 0; i < geoData_.size(); i++)
	{
		sum += geoData_[i];
	}

	return label(sum + 0.1);
}


void Foam::dataItem::show() const
{
	Info<< "Case name: " << caseName_<< nl
		<< "Case type: "<< caseType_ << nl
		<< "Case path: "<< caseDir_ << nl
		<< "Sample in x,y,z: " << m_ << " " << n_ << " " << k_ << nl
		<< "Bound box: " << xmin_ << " " << ymin_ << " " << zmin_ << " "
		<< xmax_ << " " << ymax_ << " " << zmax_ << nl
		<< "Identity: " << identity_ << nl
		<< "Ref U: " << Uref_ << nl
		<< "Ref rho: " << rhoRef_ << endl;
}


Foam::scalar Foam::dataItem::gdataSum() const
{
	scalar sum = 0;
	for (unsigned i = 0; i < geoData_.size(); i++)
	{
		sum += geoData_[i];
	}

	return sum;
}


void Foam::dataItem::setData(const word& data)
{
	wordList vtmp;

	split(vtmp, data, ',');

	if (vtmp.size() < 2)
	{
		return;
	}

    geoData_.resize(m_*n_*k_, 0);

    for (label ik = 1; ik < vtmp.size(); ik++)
	{
		wordList vtmp1;

        split(vtmp1, vtmp[ik], ':');

        label i = atoi(vtmp1[0].c_str());
        label j = atoi(vtmp1[1].c_str());
        label k = atoi(vtmp1[2].c_str());

        label ijk = i*n_*k_ + j*k_ + k;

        geoData_[ijk] = 1;
	}
}


void Foam::dataItem::setDataInfo(const word& dat)
{
	wordList vtmp;
	helyxMap::stringTokenize(vtmp, dat, ",");

	if (vtmp.size() != 12)
	{
		FatalErrorInFunction
			<< "Error set data item info, incorrect data format."
			<< exit(FatalError);
	}

	caseType_ = vtmp[0];
	caseDir_ = vtmp[1];
	caseName_ = vtmp[2];

	m_ = atoi(vtmp[3].c_str());
	n_ = atoi(vtmp[4].c_str());
	k_ = atoi(vtmp[5].c_str());

	xmin_ = atof(vtmp[6].c_str());
	ymin_ = atof(vtmp[7].c_str());
	zmin_ = atof(vtmp[8].c_str());
	xmax_ = atof(vtmp[9].c_str());
	ymax_ = atof(vtmp[10].c_str());
	zmax_ = atof(vtmp[11].c_str());
}


Foam::dataItem::dataItem(const word& dat)
{
	wordList vtmp;
	helyxMap::stringTokenize(vtmp, dat, "=");

	caseType_ = "2dcar";
	identity_ = 0;

	if (vtmp.size() == 2)
	{
		caseName_ = vtmp[0];

		wordList vtmp1;
		helyxMap::stringTokenize(vtmp1, vtmp[1], " ");

		rhoRef_ = atof(vtmp1[0].c_str());

		scalar ux = atof(vtmp1[1].c_str());
		scalar uy = atof(vtmp1[2].c_str());
		scalar uz = atof(vtmp1[3].c_str());

		Uref_[0] = ux;
		Uref_[1] = uy;
		Uref_[2] = uz;

		xmin_ = atof(vtmp1[4].c_str());
		ymin_ = atof(vtmp1[5].c_str());
		zmin_ = atof(vtmp1[6].c_str());
		xmax_ = atof(vtmp1[7].c_str());
		ymax_ = atof(vtmp1[8].c_str());
		zmax_ = atof(vtmp1[9].c_str());
	}
}


void Foam::dataItem::setData
(
	label m,
	label n,
	label k,
	const std::vector<scalar>& data
)
{
	m_ = m;
	n_ = n;
	k_ = k;
	geoData_ = data;
}


// ************************************************************************* //
