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

#include "tools/mapScalar/mapScalar.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mapScalar::mapScalar()
:
    cylindricalCoords_(false),
    f4stFormat_(false),
    mapperPtr_(nullptr),
    samplePoints_(nullptr),
    sampleTimes_(0),
    timeValue_(-1),
    startSampleTime_(-1),
    numNbrs_(5)
{}


Foam::mapScalar::mapScalar(const word& fieldName)
:
    fieldTableName_(fieldName),
    timeValue_(-1),
    numNbrs_(5)
{
    mapScalar();
}


Foam::mapScalar::mapScalar(const mapScalar& mapsc, const word& fieldName)
:
    fieldTableName_(fieldName),
    cylindricalCoords_(mapsc.cylindricalCoords_),
    f4stFormat_(mapsc.f4stFormat_),
    mapperPtr_(nullptr),
    samplePoints_(nullptr),
    sampleTimes_(0),
    timeValue_(-1),
    startSampleTime_(-1),
    numNbrs_(mapsc.numNbrs_)
{}


Foam::mapScalar::mapScalar
(
    const word& fieldName,
    const dictionary& dict,
    const scalar timeValue,
    const label neighbours
)
:
    fieldTableName_(fieldName),
    cylindricalCoords_(dict.lookupOrDefault("cylindricalCoords", false)),
    f4stFormat_(dict.lookupOrDefault("f4st", false)),
    mapperPtr_(nullptr),
    samplePoints_(nullptr),
    sampleTimes_(0),
    timeValue_(timeValue),
    startSampleTime_(-1),
    numNbrs_(neighbours)
{
    dict.readIfPresent("fieldTable", fieldTableName_);

    if (f4stFormat_)
    {
        word name_ = "surfaceData";

        char dataFile[80];
        sprintf(dataFile, "f4st/%s.f4", name_.c_str());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mapScalar::checkTable
(
    const Time& time,
    const fvPatch& patch,
    const coordinateFrame* coordFrPtr,
    bool isDefined
)
{
    // Initialise
    if (mapperPtr_.empty())
    {
        // Read cloud of points to interpolate from
        fileName samplePointsFile
        (
            time.path()
           /time.caseConstant()
           /"boundaryData"
           /patch.name()
           /"points"
        );

        samplePoints_.reset(new pointField(IFstream(samplePointsFile)()));

        if (cylindricalCoords_)
        {
            // x is r(m), y is theta(radians) and z is z(m)
            forAll(samplePoints_(), pI)
            {
                point& pntI = samplePoints_()[pI];

                scalar x = pntI.x()*Foam::cos(pntI.y());
                scalar y = pntI.x()*Foam::sin(pntI.y());

                pntI.x() = x;
                pntI.y() = y;
            }
        }

        // Allocate the interpolator
        mapperPtr_.reset
        (
            new knnInterpolation
            (
                numNbrs_,
                true,
                makePointsGlobal(samplePoints_(), coordFrPtr, isDefined)
            )
        );

        // Read the times for which data is available
        const fileName samplePointsDir = samplePointsFile.path();
        sampleTimes_ = Time::findTimes(samplePointsDir);
    }

    if (coordFrPtr && isDefined && coordFrPtr->anyDynamic())
    {
        mapperPtr_.reset
        (
            new knnInterpolation
            (
                numNbrs_,
                true,
                makePointsGlobal(samplePoints_(), coordFrPtr, isDefined)
            )
        );
    }

    // Find startSampleTime_
    label lo = -1;
    label hi = -1;
    bool foundTime = mapperPtr_().findTime
    (
        sampleTimes_,
        startSampleTime_,
        timeValue_,
        lo,
        hi
    );

    if (!foundTime)
    {
        FatalErrorInFunction
            << "Cannot find starting sampling values for current time "
            << timeValue_ << nl
            << "Have sampling values for times in directory "
            <<  time.caseConstant()/"boundaryData"/patch.name()
            << "\n    on patch " << patch.name()
            << " of field " << fieldTableName_
            << exit(FatalError);
    }

    if (lo != startSampleTime_)
    {
        startSampleTime_ = lo;
    }
}


Foam::tmp<Foam::scalarField> Foam::mapScalar::interpolateSampledData
(
    const Time& time,
    const fvPatch& patch
)
{
    // Read values
    fileName valsFile
    (
        time.path()
       /time.caseConstant()
       /"boundaryData"
       /patch.name()
       /sampleTimes_[startSampleTime_].name()
       /fieldTableName_
    );

    scalarField vals;
    IFstream(valsFile).operator()() >> vals;

    // Interpolate
    tmp<scalarField> tfld
    (
        mapperPtr_().calculateInterpolatedScalar(vals, patch.Cf())
    );

    return tfld;
}


Foam::tmp<Foam::pointField> Foam::mapScalar::makePointsGlobal
(
    pointField& samplePoints,
    const coordinateFrame* coordFrPtr,
    bool isDefined
)
{
    tmp<pointField> globalPoints(samplePoints);

    if (coordFrPtr && isDefined)
    {
        if (cylindricalCoords_)
        {
            FatalErrorInFunction
                << "Cylindrical coordinates not allowed "
                << "for use with reference frame."
                << exit(FatalError);
        }

        globalPoints = coordFrPtr->coorSys().globalPosition(samplePoints);
    }

    return globalPoints;
}


// ************************************************************************* //
