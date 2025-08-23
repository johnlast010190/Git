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

Class
    Foam::trainGeoModel

Description
    Train a machine-learning model to recognize the closest flow topology
    in the database to the one inputted.

SourceFiles
    trainGeoModel.C

\*---------------------------------------------------------------------------*/

#ifndef trainGeoModel_H
#define trainGeoModel_H

#include "sample/helyxSample/helyxSample.H"
#include "containers/dataItem/dataItem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class trainGeoModel Declaration
\*---------------------------------------------------------------------------*/

class trainGeoModel
:
    public helyxSample
{
public:

	// Public Data

		//- Enumeriate case name with number for training
		std::map<word, scalar> caseLabels_;

		std::map<word, dataItem> caseItems_;

		//- Map between search-returned ID and the ID in vCaseItems_
		std::map<label, label> idmap_;

		//- Vector version of the data items
		std::vector<dataItem> vCaseItems_;

		std::vector<word> caseNames_;

		std::vector<scalar> geoDists_;

		autoPtr<KNN> model_;


	// Constructors

        //- Default construct
        trainGeoModel();

		//- Construct from Time and dictionary
        trainGeoModel(const Time* runTime, const dictionary& dict);

		//- Construct from Time and fvMesh
		trainGeoModel(const Time* runTime, const fvMesh* mesh);


    //- Destructor
    virtual ~trainGeoModel();


	// Member Functions

		//- Clear mask array once database search tree is build
		void clearGeoData();

		void setOptions
		(
			const DynamicList<word>& argNames,
			const DynamicList<word>& argValues
		);

		void readDatabase1(const word& dbsfile);

		void readDatabase(const word& dbsfile);

		void readCaseList();

		void getCaseLabels();

		void readGeometryFile();

		void readVector(std::vector<scalar>& vect, std::ifstream& is);

		void buildModel();

		//- Build search model for specific sample resolution
		void buildModel(label m, label n, label k);

		void saveModel();

		void getSolidLabels(gridField& grid, std::vector<label>& inSolids);

		void loadModel();

		void getPrediction();

		scalar getDist(label i, label j);

		//- Get sampling data for geometry topology only
		void getGeomSample();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
