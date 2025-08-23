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
    (c) 2024-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "knn/knnInterpolation.H"
#include "src/knn.h"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::point Foam::knnInterpolation::interpolateVector
(
    const std::vector<std::vector<scalar>>& sourceV,
    const std::vector<label>& nb,
    const std::vector<scalar>& dist
)
{
    std::vector<point> prop(numNbrs_);

    for (label i = 0; i < numNbrs_; i++)
    {
        label ib = nb[i];
        std::vector<scalar> jc = sourceV[ib];

        point p;
        p.x() = jc[0];
        p.y() = jc[1];
        p.z() = jc[2];

        prop[i] = p;
    }

    scalar denom = 0.0, eps = VSMALL;
    point nom(0, 0, 0);

    for (label i = 0; i < numNbrs_; i++)
    {
        if (useInvDistSq_)
        {
            nom += prop[i]/(dist[i]*dist[i] + eps);
            denom += 1.0/(dist[i]*dist[i] + eps);
        }
        else
        {
            nom += prop[i]/(dist[i] + eps);
            denom += 1.0/(dist[i] + eps);
        }
    }

    return nom/denom;
}


Foam::point Foam::knnInterpolation::interpolateVector
(
    const std::vector<std::vector<scalar>>& sourceV,
    const scalarField& sourceAreas,
    const std::vector<label>& nb,
    const std::vector<scalar>& dist
)
{
    std::vector<point> prop(numNbrs_);
    std::vector<scalar> areas(numNbrs_);

    for (label i = 0; i < numNbrs_; i++)
    {
        label ib = nb[i];
        std::vector<scalar> jc = sourceV[ib];

        point p;
        p.x() = jc[0];
        p.y() = jc[1];
        p.z() = jc[2];

        prop[i] = p;

        // Area acts as weight in the nominator
        areas[i] = sourceAreas[ib];
    }

    scalar denom = 0.0, eps = VSMALL;
    point nom(0, 0, 0);

    for (label i = 0; i < numNbrs_; i++)
    {
        if (useInvDistSq_)
        {
            nom += (prop[i])/(dist[i]*dist[i]*areas[i]*areas[i] + eps);
            denom += (1.0)/(dist[i]*dist[i]*areas[i]*areas[i] + eps);
        }
        else
        {
            nom += (prop[i])/(dist[i]*areas[i] + eps);
            denom += (1.0)/(dist[i]*areas[i] + eps);
        }
    }

    return nom/denom;
}


Foam::point Foam::knnInterpolation::interpolateVector
(
    const std::vector<std::vector<scalar>>& sourceV,
    const pointField& sourceNormal,
    const point& ptNormal,
    const std::vector<label>& nb,
    const std::vector<scalar>& dist
)
{
    std::vector<point> prop(numNbrs_);
    std::vector<point> norm(numNbrs_);

    std::vector<scalar> dotProduct(numNbrs_);

    for (label i = 0; i < numNbrs_; i++)
    {
        label ib = nb[i];
        std::vector<scalar> jc = sourceV[ib];

        point p;
        p.x() = jc[0];
        p.y() = jc[1];
        p.z() = jc[2];

        prop[i] = p;

        // Normal vector acts as weight in the nominator
        point pnorm(sourceNormal[ib]);
        dotProduct[i] = dot(pnorm, ptNormal);
    }

    scalar denom = 0.0, eps = VSMALL;
    point nom(0, 0, 0);

    for (label i = 0; i < numNbrs_; i++)
    {
        if (useInvDistSq_)
        {
            nom +=
                (prop[i]*dotProduct[i]*dotProduct[i])/(dist[i]*dist[i] + eps);
            denom += (dotProduct[i]*dotProduct[i])/(dist[i]*dist[i] + eps);
        }
        else
        {
            nom += (prop[i]*dotProduct[i])/(dist[i] + eps);
            denom += (dotProduct[i])/(dist[i] + eps);
        }
    }

    return nom/denom;
}


Foam::scalar Foam::knnInterpolation::interpolateScalar
(
    const std::vector<scalar>& sourceV,
    const std::vector<label>& nb,
    const std::vector<scalar>& dist
)
{
    std::vector<scalar> prop(numNbrs_);

    for (label i = 0; i < numNbrs_; i++)
    {
        label ib = nb[i];
        scalar jc = sourceV[ib];
        prop[i] = jc;
    }

    scalar denom = 0.0, eps = VSMALL;
    scalar nom = 0;

    for (label i = 0; i < numNbrs_; i++)
    {
        if (useInvDistSq_)
        {
            nom += prop[i]/(dist[i]*dist[i] + eps);
            denom += 1.0/(dist[i]*dist[i] + eps);
        }
        else
        {
            nom += prop[i]/(dist[i] + eps);
            denom += 1.0/(dist[i] + eps);
        }
    }

    return nom/denom;
}


Foam::scalar Foam::knnInterpolation::calculateDistance
(
    const pointField& actualPnt,
    const triSurfaceMesh& targetSM
)
{
    scalar totalDist(0.0);

    const vector span(0.1, 0.1, 0.1);

    label found(0);

    // Loop over all actualPnt to find the nearest triangular surface
    // from targetSM
    forAll(actualPnt, pI)
    {
        const point& p(actualPnt[pI]);

        pointIndexHit inter = targetSM.nearest(p, span);
        if (inter.hit())
        {
            scalar dist(mag(inter.hitPoint() - p));
            totalDist += dist;

            found++;
        }
    }

    if (found < actualPnt.size())
    {
        Info<< "Found nearest distance for " << found
             <<"/" << actualPnt.size() << " points." << endl;
    }

    return  totalDist;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::knnInterpolation::knnInterpolation
(
    const label numNeighbours,
    bool useInverseDistanceSquared,
    const pointField& inputPoints
)
:
    numNbrs_(numNeighbours),
    useInvDistSq_(useInverseDistanceSquared),
    sourcePoints_(nullptr)
{
    sourcePoints_.set(new pointField(inputPoints));
}


Foam::knnInterpolation::knnInterpolation
(
    const label numNeighbours,
    bool useInverseDistanceSquared,
    const fileName& pointsFile
)
:
    numNbrs_(numNeighbours),
    useInvDistSq_(useInverseDistanceSquared),
    sourcePoints_(nullptr)
{
    sourcePoints_.set(new pointField(readFile(pointsFile)));
}


Foam::knnInterpolation::knnInterpolation
(
    const label numNeighbours,
    bool useInverseDistanceSquared
)
:
    numNbrs_(numNeighbours),
    useInvDistSq_(useInverseDistanceSquared),
    sourcePoints_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::knnInterpolation::calculateInterpolatedScalar
(
    const scalarField& sourceField,
    const vectorField& faceCentres
)
{
    std::vector<std::vector<scalar>> sourceP;
    field2Vector<pointField>(sourcePoints_(), sourceP);

    std::vector<scalar> sourceS;
    field2Scalar<scalarField>(sourceField, sourceS);

    KNN sourceKNN(numNbrs_);
    sourceKNN.build(sourceP);

    tmp<scalarField> tfld(new scalarField(faceCentres.size()));
    scalarField& fld = tfld.ref();

    forAll(faceCentres, cI)
    {
        const point& faceCtr = faceCentres[cI];

        std::vector<scalar> pts(3);
        pts[0] = faceCtr.x();
        pts[1] = faceCtr.y();
        pts[2] = faceCtr.z();

        std::vector<label> nb(numNbrs_);
        std::vector<scalar> dist(numNbrs_);

        sourceKNN.search(pts, nb, dist);

        scalar sc = interpolateScalar(sourceS, nb, dist);
        fld[cI] = sc;
    }

    return tfld;
}


Foam::tmp<Foam::vectorField> Foam::knnInterpolation::calculateInterpolatedVector
(
    const vectorField& sourceField,
    const vectorField& faceCentres
)
{
    std::vector<std::vector<scalar>> sourceP;
    field2Vector<pointField>(sourcePoints_(), sourceP);

    std::vector<std::vector<scalar>> sourceV;
    field2Vector<vectorField>(sourceField, sourceV);

    KNN sourceKNN(numNbrs_);
    sourceKNN.build(sourceP);

    tmp<vectorField> tfld(new vectorField(faceCentres.size()));
    vectorField& fld = tfld.ref();

    forAll(faceCentres, cI)
    {
        const point& faceCtr = faceCentres[cI];

        std::vector<scalar> pts(3);
        pts[0] = faceCtr.x();
        pts[1] = faceCtr.y();
        pts[2] = faceCtr.z();

        std::vector<label> nb(numNbrs_);
        std::vector<scalar> dist(numNbrs_);

        sourceKNN.search(pts, nb, dist);

        point p = interpolateVector(sourceV, nb, dist);
        fld[cI] = p;
    }

    return tfld;
}


Foam::tmp<Foam::pointField> Foam::knnInterpolation::readFile(const fileName& fn)
{
    tmp<pointField> ptf;

    IFstream is(fn);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << fn << nl
            << exit(FatalError);
    }

    point pt;
    while (!is.eof())
    {
        string line;
        is.getLine(line);
        if (line.length()<3)
        {
            break;
        }

        IStringStream gsp(line);
        gsp >> pt.x() >> pt.y() >> pt.z();

        ptf.ref().append(pt);
    }

    return ptf;
}


Foam::wordList Foam::knnInterpolation::timeNames(const instantList& times)
{
    wordList names(times.size());

    forAll(times, i)
    {
        names[i] = times[i].name();
    }

    return names;
}


bool Foam::knnInterpolation::findTime
(
    const instantList& times,
    const label startSampleTime,
    const scalar timeVal,
    label& lo,
    label& hi
)
{
    lo = startSampleTime;
    hi = -1;

    for (label i = startSampleTime + 1; i < times.size(); i++)
    {
        if (times[i].value() > timeVal)
        {
            break;
        }
        else
        {
            lo = i;
        }
    }

    if (lo == -1)
    {
        FatalErrorInFunction
            << "Cannot find starting sampling values for current time "
            << timeVal << nl
            << "Have sampling values for times "
            << timeNames(times) << nl
            << exit(FatalError);

        return false;
    }

    if (lo < times.size() - 1)
    {
        hi = lo + 1;
    }

    return true;
}


// ************************************************************************* //
