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

#ifndef KNN_H
#define KNN_H

// This taken verbatim from ANN.h
#ifdef WIN32
  //----------------------------------------------------------------------
  // For Microsoft Visual C++, externally accessible symbols must be
  // explicitly indicated with DLL_API, which is somewhat like "extern."
  //
  // The following ifdef block is the standard way of creating macros
  // which make exporting from a DLL simpler. All files within this DLL
  // are compiled with the DLL_EXPORTS preprocessor symbol defined on the
  // command line. In contrast, projects that use (or import) the DLL
  // objects do not define the DLL_EXPORTS symbol. This way any other
  // project whose source files include this file see DLL_API functions as
  // being imported from a DLL, wheras this DLL sees symbols defined with
  // this macro as being exported.
  //----------------------------------------------------------------------
  // TODO: Re-enable this (probably using cmake's generate_export_header()
  //       macro once we start actually doing symbol hiding!
//  #ifdef DLL_EXPORTS
//	 #define DLL_API __declspec(dllexport)
//  #else
//	#define DLL_API __declspec(dllimport)
//  #endif
  #define DLL_API
  //----------------------------------------------------------------------
  // DLL_API is ignored for all other systems
  //----------------------------------------------------------------------
#else
  #define DLL_API
#endif


#include "ANN.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

//typedef float scalar;
//typedef int label;

typedef std::vector<std::vector<double>> farray2d;
typedef std::vector<std::vector<float>> farray2s;

class DLL_API KNN
{
	public:
		int dim ;
		int m;
		int k;
		double eps ;
		int maxPts;
		bool active;
		int features; //number of features in the train vector
		int	nPts; // actual number of data points
		ANNpointArray dataPts;	// data points
		ANNpoint queryPt; // query point
		ANNidxArray	nnIdx;	// near neighbor indices
		ANNdistArray dists;	// near neighbor distances
		ANNkd_tree* kdTree;	// search structure

		KNN
		(
			const farray2d &trainData,
			int k0=2,
			double eps0=1.0e-8
		);

		inline void search
		(
			std::vector<double> &res,
			const std::vector<double> &queryVec,
			std::vector<int> &neighbour, //neighbour index in the data matrix dataPts
			std::vector<double> &vdists,
			bool interp=false
		); //search with returned results

		inline void search
		(
			std::vector<float> &res,
			const std::vector<float> &queryVec,
			std::vector<int> &neighbour, //neighbour index in the data matrix dataPts
			std::vector<float> &vdists,
			bool interp=false
		); //search with returned results

		void deactivate()
		{
			active = false;
		}

		KNN
		(
			const std::string &mdlfile
		);

		inline void build(const farray2d &trainData);//build search tree
		inline void build(const farray2s &trainData);

		void saveModel(const std::string &mdlfile);
		void loadModel(const std::string &mdlfile);

        KNN(int numNbrs=2)
        :
        k(numNbrs),
        eps(1.0e-7),
        active(true),
        dataPts(nullptr),
        queryPt(nullptr),
        nnIdx(nullptr),
        dists(nullptr),
        kdTree(nullptr)
        {
        }

		inline void search
		(
			const std::vector<double> &queryVec,
			std::vector<int> &neighbour, //neighbour index in the data matrix dataPts
			std::vector<double> &vdists
		);

        inline void search
		(
			const std::vector<double> &queryVec,
			std::vector<long int> &neighbour, //neighbour index in the data matrix dataPts
			std::vector<double> &vdists
		);


		inline void search
		(
			const std::vector<float> &queryVec,
			std::vector<int> &neighbour, //neighbour index in the data matrix dataPts
			std::vector<float> &vdists
		);


		inline void search
		(
			const std::vector<float> &queryVec,
			std::vector<long int> &neighbour, //neighbour index in the data matrix dataPts
			std::vector<float> &vdists
		);

        inline void deAllocate()
        {
            if (queryPt != NULL) annDeallocPt(queryPt);
            if (dataPts != NULL) annDeallocPts(dataPts);

            if (nnIdx != NULL)
            {
                delete [] nnIdx;
                nnIdx = NULL;
            }

            if (dists != NULL)
            {
                delete [] dists;
                dists = NULL;
            }

            if (kdTree != NULL)
            {
                delete kdTree;
                kdTree = NULL;
            }
        }

		~KNN()
        {
            deAllocate();
            annClose();
        }
};

void KNN::build
(
	const farray2d &trainData
)
{
	if (trainData.size()==0)
	{
        active = false;
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

void KNN::build
(
	const farray2s &trainData
)
{
	if (trainData.size()==0)
	{
        active = false;
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


void KNN::search
(
	const std::vector<float> &queryVec,
	std::vector<int> &neighbour, //neighbour index in the data matrix dataPts
	std::vector<float> &vdists
)
{
    if (!active) return;
    if (queryPt == NULL) return;

    for (unsigned i=0; i<queryVec.size();i++)
	{
		queryPt[i]=queryVec[i];
	}
	kdTree->annkSearch
	(
		// search
		queryPt, // query point
		k,	// number of near neighbors
		nnIdx,	// nearest neighbors (returned)
		dists,	// distance (returned)
		eps     // error bound
	);

	//put results into vector

	neighbour.resize(k,0);
	vdists.resize(k,0);

	for (int i = 0; i < k; i++)
	{
		vdists[i] = sqrtf(dists[i]); // unsquare distance
		neighbour[i]=nnIdx[i];
	}

	return;
}


void KNN::search
(
	const std::vector<float> &queryVec,
	std::vector<long int> &neighbour, //neighbour index in the data matrix dataPts
	std::vector<float> &vdists
)
{
    if (!active) return;
    if (queryPt == NULL) return;

    for (unsigned i=0; i<queryVec.size();i++)
	{
		queryPt[i]=queryVec[i];
	}
	kdTree->annkSearch
	(
		// search
		queryPt, // query point
		k,	// number of near neighbors
		nnIdx,	// nearest neighbors (returned)
		dists,	// distance (returned)
		eps     // error bound
	);

	//put results into vector

	neighbour.resize(k,0);
	vdists.resize(k,0);

	for (int i = 0; i < k; i++)
	{
		vdists[i] = sqrtf(dists[i]); // unsquare distance
		neighbour[i]=nnIdx[i];
	}

	return;
}


void KNN::search
(
	const std::vector<double> &queryVec,
	std::vector<int> &neighbour, //neighbour index in the data matrix dataPts
	std::vector<double> &vdists
)
{
    if (!active) return;
    if (queryPt == NULL) return;

    for (unsigned i=0; i<queryVec.size();i++)
	{
		queryPt[i]=queryVec[i];
	}
	kdTree->annkSearch
	(
		// search
		queryPt, // query point
		k,	// number of near neighbors
		nnIdx,	// nearest neighbors (returned)
		dists,	// distance (returned)
		eps     // error bound
	);

	//put results into vector

	neighbour.resize(k,0);
	vdists.resize(k,0);

	for (int i = 0; i < k; i++)
	{
		vdists[i] = sqrtf(dists[i]); // unsquare distance
		neighbour[i]=nnIdx[i];
	}

	return;
}

void KNN::search
(
	const std::vector<double> &queryVec,
	std::vector<long int> &neighbour, //neighbour index in the data matrix dataPts
	std::vector<double> &vdists
)
{
    if (!active) return;
    if (queryPt == NULL) return;

	for (unsigned i=0; i<queryVec.size();i++)
	{
		queryPt[i]=queryVec[i];
	}
	kdTree->annkSearch
	(
		// search
		queryPt, // query point
		k,	// number of near neighbors
		nnIdx,	// nearest neighbors (returned)
		dists,	// distance (returned)
		eps     // error bound
	);

	//put results into vector

	neighbour.resize(k,0);
	vdists.resize(k,0);

	for (int i = 0; i < k; i++)
	{
		vdists[i] = sqrtf(dists[i]); // unsquare distance
		neighbour[i]=nnIdx[i];
	}

	return;
}


void KNN::search
(
	std::vector<double> &res,
	const std::vector<double> &queryVec,
	std::vector<int> &nb,
	std::vector<double> &vdists,
	bool interp
)
{
    if (!active) return;

	res.resize(dim);
	search(queryVec,nb,vdists);
	if (!interp)
	{
		int ic=nb[0];
		for (int j=0;j<dim;j++)
		{
			res[j]=dataPts[ic][j];
		}
	}
	else
	{
		double small=1.0e-30;
		for (int j=0;j<dim;j++)
		{
			double sum=0;
			double dsum=0;
			for (int ik=0;ik<k;ik++)
			{
				int ick=nb[ik];
				double ditem=1.0/(small+vdists[ik]);
				sum+=dataPts[ick][j]*ditem;
				dsum+=ditem;
			}
			res[j]=sum/dsum;
		}
	}

}

void KNN::search
(
	std::vector<float> &res,
	const std::vector<float> &queryVec,
	std::vector<int> &nb,
	std::vector<float> &vdists,
	bool interp
)
{
	res.resize(dim);
	search(queryVec,nb,vdists);
	if (!interp)
	{
		int ic=nb[0];
		for (int j=0;j<dim;j++)
		{
			res[j]=dataPts[ic][j];
		}
	}
	else
	{
		double small=1.0e-30;
		for (int j=0;j<dim;j++)
		{
			double sum=0;
			double dsum=0;
			for (int ik=0;ik<k;ik++)
			{
				int ick=nb[ik];
				double ditem=1.0/(small+vdists[ik]);
				sum+=dataPts[ick][j]*ditem;
				dsum+=ditem;
			}
			res[j]=sum/dsum;
		}
	}

}

#endif // KNN_H
