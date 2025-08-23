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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2017-2022 Engys Ltd.

Description
    Calulate the face centres and areas.

    Calculate the centre by breaking the face into triangles using the face
    centre and area-weighted averaging their centres.  This method copes with
    small face-concavity.

\*---------------------------------------------------------------------------*/

#include "meshes/primitiveMesh/primitiveMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcFaceCentresAndAreas() const
{
    if (debug)
    {
        Pout<< "primitiveMesh::calcFaceCentresAndAreas() : "
            << "Calculating face centres and face areas"
            << endl;
    }

    // It is an error to attempt to recalculate faceCentres
    // if the pointer is already set
    if (faceCentresPtr_ || faceAreasPtr_ || magFaceAreasPtr_)
    {
        FatalErrorInFunction
            << "Face centres or face areas already calculated"
            << abort(FatalError);
    }

    faceCentresPtr_ = new vectorField(nFaces());
    vectorField& fCtrs = *faceCentresPtr_;

    faceAreasPtr_ = new vectorField(nFaces());
    vectorField& fAreas = *faceAreasPtr_;

    magFaceAreasPtr_ = new scalarField(nFaces());
    scalarField& magfAreas = *magFaceAreasPtr_;

    makeFaceCentresAndAreas(points(), fCtrs, fAreas, magfAreas);

    if (debug)
    {
        Pout<< "primitiveMesh::calcFaceCentresAndAreas() : "
            << "Finished calculating face centres and face areas"
            << endl;
    }
}


void Foam::primitiveMesh::faceAreaAndCentre
(
    const face& f,
    const pointField& p,
    vector& fCentre,
    vector& fArea
)
{
 #if defined(HELYX_SP)
    Vector<double> fC = Zero;
    Vector<double> fA = Zero;

    label nPoints = f.size();
    Vector<double> pt0 = p[f[0]];
    Field<Vector<double>> facePts(nPoints,Zero);
    // Translate by face by point 0
    for (label fp = 1; fp < nPoints; fp++)
    {
       facePts[fp] = p[f[fp]];
       facePts[fp] -= pt0;
    }

    // If the face is a triangle, do a direct calculation for efficiency
    // and to avoid round-off error-related problems
    if (nPoints == 3)
    {
        fC = (1.0/3.0)*(facePts[0] + facePts[1] + facePts[2]);
        fA =
            0.5*((facePts[1] - facePts[0])^(facePts[2] - facePts[0]));
    }
    else
    {
        Vector<double> sumN = Zero;
        double sumA = 0.0;
        Vector<double> sumAc = Zero;

        Vector<double> midPt = Zero;
        forAll(f, pI)
        {
            midPt += facePts[pI];
        }
        midPt /= nPoints;

        Field<Vector<double>> n(nPoints);
        forAll(f, pI)
        {
            const Vector<double>& pt = facePts[pI];
            const Vector<double>& pNext = facePts[f.fcIndex(pI)];
            n[pI] = (pNext - pt)^(midPt - pt);
            sumN += n[pI];
        }

        const Vector<double> sumNHat = normalised(sumN);
        forAll(f, pI)
        {
            const Vector<double>& pt = facePts[pI];
            const Vector<double>& pNext = facePts[f.fcIndex(pI)];
            const Vector<double> c = pt + pNext + midPt;
            double a = n[pI] & sumNHat;
            sumA += a;
            sumAc += a*c;
        }

        // This is to deal with zero-area faces. Mark very small faces
        // to be detected in e.g., processorPolyPatch.
        if (sumA > VSMALL)
        {
            fC = (1.0/3.0)*sumAc/sumA;
        }
        else
        {
            fC = midPt;
        }
        fA = 0.5*sumN;
    }

    //Correct face centres for face point translations
    fC += pt0;
    fCentre = fC;
    fArea = fA;
#else
    fCentre = Zero;
    fArea = Zero;

    label nPoints = f.size();
    vector pt0 = p[f[0]];
    pointField facePts(nPoints,Zero);
    // Translate by face by point 0
    for (label fp = 1; fp < nPoints; fp++)
    {
       facePts[fp] = p[f[fp]];
       facePts[fp] -= pt0;
    }

    // If the face is a triangle, do a direct calculation for efficiency
    // and to avoid round-off error-related problems
    if (nPoints == 3)
    {
        fCentre = (1.0/3.0)*(facePts[0] + facePts[1] + facePts[2]);
        fArea =
            0.5*((facePts[1] - facePts[0])^(facePts[2] - facePts[0]));
    }
    else
    {
        vector sumN = Zero;
        scalar sumA = 0.0;
        vector sumAc = Zero;

        point midPt = Zero;
        forAll(f, pI)
        {
            midPt += facePts[pI];
        }
        midPt /= nPoints;

        pointField n(nPoints);
        forAll(f, pI)
        {
            const point& pt = facePts[pI];
            const point& pNext = facePts[f.fcIndex(pI)];
            n[pI] = (pNext - pt)^(midPt - pt);
            sumN += n[pI];
        }

        const vector sumNHat = normalised(sumN);
        forAll(f, pI)
        {
            const point& pt = facePts[pI];
            const point& pNext = facePts[f.fcIndex(pI)];
            const vector c = pt + pNext + midPt;
            scalar a = n[pI] & sumNHat;
            sumA += a;
            sumAc += a*c;
        }

        // This is to deal with zero-area faces. Mark very small faces
        // to be detected in e.g., processorPolyPatch.
        if (sumA > VSMALL)
        {
            fCentre = (1.0/3.0)*sumAc/sumA;
        }
        else
        {
            fCentre = midPt;
        }
        fArea = 0.5*sumN;
    }
    //Correct face centres for face point translations
    fCentre += pt0;
#endif
}


void Foam::primitiveMesh::makeFaceCentresAndAreas
(
    const pointField& p,
    vectorField& fCtrs,
    vectorField& fAreas,
    scalarField& magfAreas
) const
{
    const faceList& fs = faces();
    forAll(fs, facei)
    {
        const face& f = fs[facei];
        faceAreaAndCentre(f,p,fCtrs[facei],fAreas[facei]);
        magfAreas[facei] = max(mag(fAreas[facei]), VSMALL);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// TODO: Inline me
const Foam::scalarField& Foam::primitiveMesh::magFaceAreas() const
{
    if (!magFaceAreasPtr_)
    {
        calcFaceCentresAndAreas();
    }

    return *magFaceAreasPtr_;
}

// ************************************************************************* //
