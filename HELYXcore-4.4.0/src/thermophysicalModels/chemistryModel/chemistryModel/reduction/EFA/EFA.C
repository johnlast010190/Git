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
    (c) 2016-2021 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "chemistryModel/reduction/EFA/EFA.H"
#include "chemistryModel/reduction/EFA/SortableListEFA.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::chemistryReductionMethods::EFA<ThermoType>::EFA
(
    const IOdictionary& dict,
    chemistryModel<ThermoType>& chemistry
)
:
    chemistryReductionMethod<ThermoType>(dict, chemistry),
    sC_(this->nSpecie(),0),
    sH_(this->nSpecie(),0),
    sO_(this->nSpecie(),0),
    sN_(this->nSpecie(),0),
    sortPart_(0.05)
{
    for (label i=0; i<this->nSpecie(); i++)
    {
        const List<specieElement>& curSpecieComposition =
            chemistry.mixture().specieComposition(i);

        // for all elements in the current species
        forAll(curSpecieComposition, j)
        {
            const specieElement& curElement =
            curSpecieComposition[j];
            if (curElement.name() == "C")
            {
                sC_[i] = curElement.nAtoms();
            }
            else if (curElement.name() == "H")
            {
                sH_[i] = curElement.nAtoms();
            }
            else if (curElement.name() == "O")
            {
                sO_[i] = curElement.nAtoms();
            }
            else if (curElement.name() == "N")
            {
                sN_[i] = curElement.nAtoms();
            }
            else
            {
                Info<< "element not considered"<< endl;
            }
        }
    }
    if (this->coeffsDict_.found("sortPart"))
    {
        sortPart_ = this->coeffsDict_.template lookup<scalar>("sortPart");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::chemistryReductionMethods::EFA<ThermoType>::~EFA()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::chemistryReductionMethods::EFA<ThermoType>::reduceMechanism
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    List<label>& ctos,
    DynamicList<label>& stoc,
    const label li
)
{
    chemistryReductionMethod<ThermoType>::initReduceMechanism();

    scalarField c1(this->chemistry_.nEqns(), 0.0);

    for (label i=0; i<this->nSpecie(); i++)
    {
        c1[i] = c[i];
    }

    c1[this->nSpecie()] = T;
    c1[this->nSpecie() + 1] = p;


    // Number of initialised rAB for each lines
    Field<label> NbrABInit(this->nSpecie(),0);

    // Position of the initialised rAB, -1 when not initialised
    RectangularMatrix<label> rABPos(this->nSpecie(), this->nSpecie(), -1);
    RectangularMatrix<scalar> CFluxAB(this->nSpecie(), this->nSpecie(), 0.0);
    RectangularMatrix<scalar> HFluxAB(this->nSpecie(), this->nSpecie(), 0.0);
    RectangularMatrix<scalar> OFluxAB(this->nSpecie(), this->nSpecie(), 0.0);
    RectangularMatrix<scalar> NFluxAB(this->nSpecie(), this->nSpecie(), 0.0);
    scalar CFlux(0.0), HFlux(0.0), OFlux(0.0), NFlux(0.0);
    label nbPairs(0);

    // Index of the other species involved in the rABNum
    RectangularMatrix<label> rABOtherSpec(this->nSpecie(), this->nSpecie(), -1);

    forAll(this->chemistry_.reactions(), i)
    {
        const Reaction<ThermoType>& R = this->chemistry_.reactions()[i];
        // for each reaction compute omegai
        scalar omegaf, omegar;
        R.omega(p, T, c1, li, omegaf, omegar);

        scalar fr = mag(omegaf) + mag(omegar);
        scalar NCi(0.0), NHi(0.0), NOi(0.0), NNi(0.0);
        forAll(R.lhs(),s)
        {
            label curIndex = R.lhs()[s].index;
            scalar stoicCoeff = R.lhs()[s].stoichCoeff;
            NCi += sC_[curIndex]*stoicCoeff;
            NHi += sH_[curIndex]*stoicCoeff;
            NOi += sO_[curIndex]*stoicCoeff;
            NNi += sN_[curIndex]*stoicCoeff;
        }
        // element conservation means that total number of element is
        // twice the number on the left
        NCi *=2;
        NHi *=2;
        NOi *=2;
        NNi *=2;
        // then for each source/sink pairs, compute the element flux
        forAll(R.lhs(), Ai)// compute the element flux
        {
            label A = R.lhs()[Ai].index;
            scalar Acoeff = R.lhs()[Ai].stoichCoeff;

            forAll(R.rhs(), Bi)
            {
                label B = R.rhs()[Bi].index;
                scalar Bcoeff = R.rhs()[Bi].stoichCoeff;
                // if a pair in the reversed way has not been initialised
                if (rABPos(B, A) == -1)
                {
                    label otherS = rABPos(A, B);
                    if (otherS == -1)
                    {
                        rABPos(A, B) = NbrABInit[A]++;
                        otherS = rABPos(A, B);
                        nbPairs++;
                        rABOtherSpec(A, otherS) = B;
                    }
                    if (NCi>VSMALL)
                    {
                        CFluxAB(A, otherS) +=
                            fr*Acoeff*sC_[A]*Bcoeff*sC_[B]/NCi;
                    }
                    if (NHi>VSMALL)
                    {
                        HFluxAB(A, otherS) +=
                            fr*Acoeff*sH_[A]*Bcoeff*sH_[B]/NHi;
                    }
                    if (NOi>VSMALL)
                    {
                        OFluxAB(A, otherS) +=
                            fr*Acoeff*sO_[A]*Bcoeff*sO_[B]/NOi;
                    }
                    if (NNi>VSMALL)
                    {
                        NFluxAB(A, otherS) +=
                            fr*Acoeff*sN_[A]*Bcoeff*sN_[B]/NNi;
                    }
                }
                // If a pair BA is initialised,
                // add the element flux to this pair
                else
                {
                    label otherS = rABPos(B, A);
                    if (NCi>VSMALL)
                    {
                        CFluxAB(B, otherS) +=
                            fr*Acoeff*sC_[A]*Bcoeff*sC_[B]/NCi;
                    }
                    if (NHi>VSMALL)
                    {
                        HFluxAB(B, otherS) +=
                            fr*Acoeff*sH_[A]*Bcoeff*sH_[B]/NHi;
                    }
                    if (NOi>VSMALL)
                    {
                        OFluxAB(B, otherS) +=
                            fr*Acoeff*sO_[A]*Bcoeff*sO_[B]/NOi;
                    }
                    if (NNi>VSMALL)
                    {
                        NFluxAB(B, otherS) +=
                            fr*Acoeff*sN_[A]*Bcoeff*sN_[B]/NNi;
                    }
                }
                if (NCi>VSMALL)
                {
                    CFlux += fr*Acoeff*sC_[A]*Bcoeff*sC_[B]/NCi;
                }
                if (NHi>VSMALL)
                {
                    HFlux += fr*Acoeff*sH_[A]*Bcoeff*sH_[B]/NHi;
                }
                if (NOi>VSMALL)
                {
                    OFlux += fr*Acoeff*sO_[A]*Bcoeff*sO_[B]/NOi;
                }
                if (NNi>VSMALL)
                {
                    NFlux += fr*Acoeff*sN_[A]*Bcoeff*sN_[B]/NNi;
                }
            }
        }
    }

    // Select species according to the total flux cutoff (1-tolerance)
    // of the flux is included
    for (label i=0; i<this->nSpecie(); i++)
    {
        this->activeSpecies_[i] = false;
    }

    if (CFlux > VSMALL)
    {
        SortableListEFA<scalar> pairsFlux(nbPairs);
        labelList source(nbPairs);
        labelList sink(nbPairs);
        label nP(0);
        for (int A=0; A<this->nSpecie() ; A++)
        {
            for (int j=0; j<NbrABInit[A]; j++)
            {
                label pairIndex = nP++;
                pairsFlux[pairIndex] = 0.0;
                label B = rABOtherSpec(A, j);
                pairsFlux[pairIndex] += CFluxAB(A, j);
                source[pairIndex] = A;
                sink[pairIndex] = B;
            }
        }

        // Sort in descending orders the source/sink pairs until the cutoff is
        // reached. The sorting is done on part of the pairs because a small
        // part of these pairs are responsible for a large part of the element
        // flux
        scalar cumFlux(0.0);
        scalar threshold((1-this->tolerance())*CFlux);
        label startPoint(0);
        label nbToSort(static_cast<label> (nbPairs*sortPart_));
        nbToSort = max(nbToSort,1);

        bool cumRespected(false);
        while (!cumRespected)
        {
            if (startPoint >= nbPairs)
            {
                FatalErrorInFunction
                    << "startPoint outside number of pairs without reaching"
                    << "100% flux"
                    << exit(FatalError);
            }

            label nbi = min(nbToSort, nbPairs-startPoint);
            pairsFlux.partialSort(nbi, startPoint);
            const labelList& idx = pairsFlux.indices();
            for (int i=0; i<nbi; i++)
            {
                cumFlux += pairsFlux[idx[startPoint+i]];
                if (!this->activeSpecies_[source[idx[startPoint+i]]])
                {
                    this->activeSpecies_[source[idx[startPoint+i]]] = true;
                }
                if (!this->activeSpecies_[sink[idx[startPoint+i]]])
                {
                    this->activeSpecies_[sink[idx[startPoint+i]]] = true;
                }
                if (cumFlux >= threshold)
                {
                    cumRespected = true;
                    break;
                }
            }
            startPoint += nbToSort;
        }
    }

    if (HFlux > VSMALL)
    {
        SortableListEFA<scalar> pairsFlux(nbPairs);
        labelList source(nbPairs);
        labelList sink(nbPairs);
        label nP(0);
        for (int A=0; A<this->nSpecie() ; A++)
        {
            for (int j=0; j<NbrABInit[A]; j++)
            {
                label pairIndex = nP++;
                pairsFlux[pairIndex] = 0.0;
                label B = rABOtherSpec(A, j);
                pairsFlux[pairIndex] += HFluxAB(A, j);
                source[pairIndex] = A;
                sink[pairIndex] = B;
            }
        }

        scalar cumFlux(0.0);
        scalar threshold((1-this->tolerance())*HFlux);
        label startPoint(0);
        label nbToSort(static_cast<label> (nbPairs*sortPart_));
        nbToSort = max(nbToSort,1);

        bool cumRespected(false);
        while (!cumRespected)
        {
            if (startPoint >= nbPairs)
            {
                FatalErrorInFunction
                    << "startPoint outside number of pairs without reaching"
                    << "100% flux"
                    << exit(FatalError);
            }

            label nbi = min(nbToSort, nbPairs-startPoint);
            pairsFlux.partialSort(nbi, startPoint);
            const labelList& idx = pairsFlux.indices();
            for (int i=0; i<nbi; i++)
            {
                cumFlux += pairsFlux[idx[startPoint+i]];

                if (!this->activeSpecies_[source[idx[startPoint+i]]])
                {
                    this->activeSpecies_[source[idx[startPoint+i]]] = true;
                }
                if (!this->activeSpecies_[sink[idx[startPoint+i]]])
                {
                    this->activeSpecies_[sink[idx[startPoint+i]]] = true;
                }
                if (cumFlux >= threshold)
                {
                    cumRespected = true;
                    break;
                }
            }
            startPoint += nbToSort;
        }
    }

    if (OFlux > VSMALL)
    {
        SortableListEFA<scalar> pairsFlux(nbPairs);
        labelList source(nbPairs);
        labelList sink(nbPairs);
        label nP(0);
        for (int A=0; A<this->nSpecie() ; A++)
        {
            for (int j=0; j<NbrABInit[A]; j++)
            {
                label pairIndex = nP++;
                pairsFlux[pairIndex] = 0.0;
                label B = rABOtherSpec(A, j);
                pairsFlux[pairIndex] += OFluxAB(A, j);
                source[pairIndex] = A;
                sink[pairIndex] = B;
            }
        }
        scalar cumFlux(0.0);
        scalar threshold((1-this->tolerance())*OFlux);
        label startPoint(0);
        label nbToSort(static_cast<label> (nbPairs*sortPart_));
        nbToSort = max(nbToSort,1);

        bool cumRespected(false);
        while (!cumRespected)
        {
            if (startPoint >= nbPairs)
            {
                FatalErrorInFunction
                    << "startPoint outside number of pairs without reaching"
                    << "100% flux"
                    << exit(FatalError);
            }

            label nbi = min(nbToSort, nbPairs-startPoint);
            pairsFlux.partialSort(nbi, startPoint);
            const labelList& idx = pairsFlux.indices();
            for (int i=0; i<nbi; i++)
            {
                cumFlux += pairsFlux[idx[startPoint+i]];

                if (!this->activeSpecies_[source[idx[startPoint+i]]])
                {
                    this->activeSpecies_[source[idx[startPoint+i]]] = true;
                }
                if (!this->activeSpecies_[sink[idx[startPoint+i]]])
                {
                    this->activeSpecies_[sink[idx[startPoint+i]]] = true;
                }
                if (cumFlux >= threshold)
                {
                    cumRespected = true;
                    break;
                }
            }
            startPoint += nbToSort;
        }
    }

    if (NFlux > VSMALL)
    {
        SortableListEFA<scalar> pairsFlux(nbPairs);
        labelList source(nbPairs);
        labelList sink(nbPairs);
        label nP(0);
        for (int A=0; A<this->nSpecie() ; A++)
        {
            for (int j=0; j<NbrABInit[A]; j++)
            {
                label pairIndex = nP++;
                pairsFlux[pairIndex] = 0.0;
                label B = rABOtherSpec(A, j);
                pairsFlux[pairIndex] += NFluxAB(A, j);
                source[pairIndex] = A;
                sink[pairIndex] = B;
            }
        }
        scalar cumFlux(0.0);
        scalar threshold((1-this->tolerance())*NFlux);
        label startPoint(0);
        label nbToSort(static_cast<label> (nbPairs*sortPart_));
        nbToSort = max(nbToSort,1);

        bool cumRespected(false);
        while (!cumRespected)
        {
            if (startPoint >= nbPairs)
            {
                FatalErrorInFunction
                    << "startPoint outside number of pairs without reaching"
                    << "100% flux"
                    << exit(FatalError);
            }

            label nbi = min(nbToSort, nbPairs-startPoint);
            pairsFlux.partialSort(nbi, startPoint);
            const labelList& idx = pairsFlux.indices();
            for (int i=0; i<nbi; i++)
            {
                cumFlux += pairsFlux[idx[startPoint+i]];

                if (!this->activeSpecies_[source[idx[startPoint+i]]])
                {
                    this->activeSpecies_[source[idx[startPoint+i]]] = true;
                }
                if (!this->activeSpecies_[sink[idx[startPoint+i]]])
                {
                    this->activeSpecies_[sink[idx[startPoint+i]]] = true;
                }
                if (cumFlux >= threshold)
                {
                    cumRespected = true;
                    break;
                }
            }
            startPoint += nbToSort;
        }
    }

    chemistryReductionMethod<ThermoType>::endReduceMechanism(ctos, stoc);
}


// ************************************************************************* //
