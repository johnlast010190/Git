/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |
|    o     o     |  ENGYS Ltd. <http://engys.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of a specialisation of the Approximate Nearest Neighbor
    Library (ANN) library for HELYXcore.  This software is provided under the
    provisions of the Lesser GNU Public License (LGPL), in accordance with the
    requirements of the ANN library.  See the file ../ReadMe.txt for further
    information.

Copyright
    (c) 2020 Engys Ltd.

\*---------------------------------------------------------------------------*/


#include "knn.h"


void KNN::saveModel(const std::string &mdlfile)
{
	//save search tree to model file
	std::ofstream os(mdlfile.c_str(),std::ios::out);
	if (!os)
	{
		std::cout
			<< "Error save KNN model into file, file\n"
			<< mdlfile<<" cannot open for writing."<<std::endl;
		std::exit(1);
	}
	os<<dim<<" "<<m<<" "<<k<<" "<<eps<<" "<<maxPts<<" "<<features
	<<" "<<nPts<<std::endl;
	kdTree->Dump
	(	// dump entire tree
		ANNtrue,
		os
	);
	return;
}

KNN::KNN(const std::string &mdlfile)
{
	active=true;
	std::ifstream is(mdlfile.c_str(),std::ios::in);
	if (!is)
	{
		std::cout
			<< "Error build KNN from model file, file\n"
			<< mdlfile<<" cannot open for reading."<<std::endl;
		std::exit(1);
	}
	is>>dim>>m>>k>>eps>>maxPts>>features>>nPts;
	kdTree=new ANNkd_tree(is);
}

/*
KNN::KNN()
{
	k =2;
	eps =1.0e-7;
	active=true;
}


void KNN::build
(
	const farray2d &trainData
)
{
	if (trainData.size()==0)
	{
		return;
	}

	dim=trainData[0].size();

	m =trainData.size();
	features=dim;
	maxPts=trainData.size(); //actual number of data points
	queryPt = annAllocPt(dim);					// allocate query point
	dataPts = annAllocPts(maxPts, dim);			// allocate data points
	nnIdx = new ANNidx[k];						// allocate near neigh indices
	dists = new ANNdist[k];						// allocate near neighbor dists

	nPts = 0;									// read data points
	//get data
	//Info<<"dim,m,maxpts:"<<	maxPts<<" "<<dim<<" "<<	trainData.size()<<endl;
	for (int i=0;i<maxPts;i++)
	{
		for (int j=0; j<dim;j++)
		{
			dataPts[i][j]=trainData[i][j];
		}
	}

	int nPts=trainData.size();
	nPts=maxPts;

    kdTree = new ANNkd_tree
    (
       // build search structure
		dataPts,	// the data points
		nPts, // number of points
		dim
	); // dimension of space


}

*/

KNN::KNN
(
	const farray2d &trainData,
	int k0,
	double eps0
)
{
	k =k0;
	active=true;
	eps =eps0;
	build(trainData);
}

/*
KNN::~KNN()
{

//	delete [] nnIdx; // clean things up
//    delete [] dists;
 //   delete kdTree;
	annClose();	// done with ANN
}
*/


