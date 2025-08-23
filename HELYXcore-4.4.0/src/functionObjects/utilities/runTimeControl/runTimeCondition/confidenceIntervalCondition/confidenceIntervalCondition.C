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
    (c) 2018, 2020 Engys Ltd
    (c) 2015 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/Fields/Field/SubField.H"
#include "runTimeControl/runTimeCondition/confidenceIntervalCondition/confidenceIntervalCondition.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/etcFiles/etcFiles.H"
#include "primitives/functions/Function1/Table/TableBase.H"
#include "primitives/functions/Function1/Table/Table.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "matrices/simpleMatrix/simpleMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeControls
{
    defineTypeNameAndDebug(confidenceIntervalCondition, 0);
    addToRunTimeSelectionTable
    (
        runTimeCondition,
        confidenceIntervalCondition,
        dictionary
    );
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::runTimeControls
    ::confidenceIntervalCondition::intvTypes,
    2
>::names[] =
{
    "relative",
    "absolute"
};

const Foam::NamedEnum
<
    Foam::functionObjects::runTimeControls
    ::confidenceIntervalCondition::intvTypes,
    2
> Foam::functionObjects::runTimeControls
    ::confidenceIntervalCondition::intvTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::runTimeControls::confidenceIntervalCondition
::calc
(
    bool& satisfied
)
{
    scalar equilibratedAverage = 0;
    scalar equilibratedUncertainty = 1e-5;
    scalar tCrit = 1;
    scalar transientTime = 0;

    List<scalar> bufferX(X_.size()-filterRadius_);
    List<scalar> bufferT(T_.size()-filterRadius_);
    forAll(bufferX, i)
    {
        bufferX[i] = X_[i];
        bufferT[i] = T_[i];
    }

    List<scalar> smoothX(bufferX);
    List<scalar> highPass(bufferX);
    if (filterRadius_!=0 && 2*filterRadius_ < X_.size())
    {
        smoothX = lowPassFilter(X_);
    }

    label lastTransientSample = 0;
    forAll(bufferT, i)
    {
       if (bufferT[i]>minTransientTime_)
       {
            lastTransientSample = i;
            break;
       }
    }

    label cutOffSize = bufferX.size() - lastTransientSample;
    List<scalar> truncT(cutOffSize);
    scalarField truncX(cutOffSize);
    for (int i=0; i<cutOffSize; ++i)
    {
        label oldIndex = i + lastTransientSample;
        truncT[i] = bufferT[oldIndex];
        truncX[i] = bufferX[oldIndex];
    }

    transientTime = findGlobalMinimum
    (
        (smoothValueForTransientDetection_) ? smoothX : truncX,
        truncT,
        truncX.size()
    );

    transientTime = max(transientTime, minTransientTime_);

    scalar minDistance = GREAT;
    label transientLabel = 0;

    forAll(bufferT, i)
    {
        scalar distance = mag(bufferT[i]-transientTime);
        if (distance<minDistance)
        {
            minDistance = mag(distance);
            transientLabel = i;
        }
    }

    // Ensure there are at least 2 steady samples
    transientLabel = min(transientLabel, mag(bufferX.size()-2));

    label steadyStateSize = bufferX.size() - transientLabel;

    List<scalar> ssT(steadyStateSize);
    scalarField ssX(steadyStateSize);
    scalarField ssSmoothX(steadyStateSize);

    for (int i=0; i<steadyStateSize; ++i)
    {
        label oldIndex = i + transientLabel;
        ssT[i] = bufferT[oldIndex];
        ssX[i] = bufferX[oldIndex];
        ssSmoothX[i] = smoothX[oldIndex];
    }

    label sampleSize = ssX.size();

    List<scalar> aCorr = autoCorrelationFactor(ssX);

    scalar aCorrTime = autoCorrelationTime( aCorr);

    //- We don't include the first element since it is always 1 and
    //- for 0 auto-correlation Neff should equal N.
    label endAutoCorr = min(aCorrTime, aCorr.size()-1);

    scalar sumOfAutoCorr = 0;
    for (int i=1; i<=endAutoCorr; ++i)
    {
        sumOfAutoCorr += aCorr[i];
    }

    scalar nEff =
    (
        (sampleSize - 2*aCorrTime- 1)
      + (static_cast<scalar>(aCorrTime*(aCorrTime+1))/sampleSize)
    )/(1+2*sumOfAutoCorr);

    nEff += 1;
    nEff = max(nEff, 2);

    equilibratedAverage = meanData(ssX, 0, ssX.size());

    scalar variance = 0;
    forAll(ssX, i)
    {
        variance += sqr(ssX[i] - equilibratedAverage);
    }

    equilibratedUncertainty =
        sqrt(nEff/(sampleSize*(nEff-1))*variance)/sqrt(nEff);

    label DoF = std::ceil(nEff) -1;

//    computeEquilibratedUncertainty(ssX, aCorrTime, equilibratedAverage, equilibratedUncertainty, DoF);

    tCrit = tDistribution_->value(confidence_, DoF);

    scalar tolerance = interval_;

    if (intervalType_ == itRel)
    {
        tolerance *= mag(equilibratedAverage);
    }

    //Estimate Sample Size
    scalar currentStd = sqrt(variance/(nEff));
    scalar zCrit = tDistribution_->value(confidence_, 1e04);
    scalar minimumSamplesNeeded = std::ceil(sqr(zCrit*currentStd/tolerance));

    scalar confIntv = tCrit*equilibratedUncertainty;

    Log << tab << "value                    - " << smoothX.last() << nl
        << tab << "N samples                - " << bufferX.size() << nl
        << tab << "Time                     - " << bufferT.last() << nl
        << tab << "Transient Time           - " << transientTime << nl
        << tab << "Steady-State samples     - " << ssX.size() << nl
        << tab << "Samples needed (estim.)  - " << minimumSamplesNeeded << nl
        << tab << "Degrees of Freedom       - " << DoF + 1 << nl
        << tab << "Auto-Correlation Int     - " << sumOfAutoCorr<< nl
        << tab << "mean                     - " << equilibratedAverage << nl
        << tab << "std error                - " << equilibratedUncertainty << nl
        << tab << "confidence level         - "
        << confidence_ * 100.0 << "%" << nl
        << tab << "target confidence interval  - " << tolerance << nl
        << tab << "current confidence interval - " << confIntv << nl<<endl;


    scalar slope = getRegressedDataSlope
    (
        ssT,
        (smoothValueForConvergence_) ? ssSmoothX : ssX
    );

    satisfied = checkConvergence
    (
        tolerance,
        confIntv,
        slope,
        ssX.size()
    );

    if (Pstream::master() || !Pstream::parRun())
    {
        scalar upperBound = equilibratedAverage +
            (equilibratedUncertainty*tCrit);
        scalar lowerBound = equilibratedAverage -
            (equilibratedUncertainty*tCrit);
        file() << bufferX.size() << tab << bufferT.last() << tab <<
        transientTime << tab
               << bufferX.last() << tab << equilibratedAverage << tab
               << upperBound << tab << lowerBound << tab << smoothX.last()
               << tab << slope << endl;
    }
}

Foam::List<Foam::scalar> Foam::functionObjects::runTimeControls
::confidenceIntervalCondition
::lowPassFilter
(
    const List<Foam::scalar>& X
)
{

    scalarField smoothedX(X.size()-filterRadius_);

    label kernelSize = std::ceil(2*filterRadius_+1);
    scalar sigma = filterRadius_/2.;
    scalarField kernel(kernelSize, 0);
    scalar normFactor = 1.0 / (sqrt(2.0 * M_PI) * sigma);
    forAll(kernel, i)
    {
        scalar x = i - filterRadius_;
        kernel[i] = exp(-(x*x)/(2*sigma*sigma))/normFactor;
    }

    for (int i = 0; i < smoothedX.size(); i++)
    {
        scalar sum = 0.0;
        scalar weightSum = 0.0;
        for (int j = 0; j < kernelSize; j++)
        {
            int index = i + j - filterRadius_ ;
            if (index >= 0 && index < X.size())
            {
                sum += kernel[j] * X[index];
                weightSum += kernel[j];
            }
        }
        smoothedX[i] = sum/weightSum;
    }

    return std::move(smoothedX);
}

Foam::scalar Foam::functionObjects::runTimeControls::confidenceIntervalCondition
::findGlobalMinimum
(
    const Foam::List<scalar>& X,
    const Foam::List<scalar>& t,
    const Foam::label& k
)
{
    if (k<9)
    {
        return 0;
    }

    label h = k;
    scalar tmin = 0.0;

    List<scalar> Xnew = X;

    // Iterate until h becomes 2.
    List<scalar> XrevSE = reverseStandardError(Xnew, h);
    scalar meanTime = 0;

    DynamicList<scalar> minTimes;
    while (h > 8)
    {
        // Identify indices i (with 0 < i < h-1) that are local minima of XrevSE
        // and for which the spread in Xnew[i-1:i+2] exceeds max_err.
        std::vector<label> idx_gut;
        for (label i = 1; i < h - 1; ++i)
        {
            if (XrevSE[i] < XrevSE[i - 1] && XrevSE[i] < XrevSE[i + 1])
            {
                // Compute local maximum and minimum in Xnew[i-1], Xnew[i], Xnew[i+1]
                scalar local_max = std::max({Xnew[i - 1], Xnew[i], Xnew[i + 1]});
                scalar local_min = std::min({Xnew[i - 1], Xnew[i], Xnew[i + 1]});
                if (mag(local_max - local_min) > 1e-5)
                {
                    idx_gut.push_back(i);
                }
            }
        }

        if (debug && k==1000)
        {
            word name("reverseError_"+std::to_string(h));
            fileName fileN(name);
            OFstream os(fileN);
            forAll(XrevSE, i)
            {
                os << XrevSE[i] << nl;
            }
        }


        label idx = 0;
        if (!idx_gut.empty())
        {
            idx = idx_gut[0];
            {
                // Choose the index among idx_gut that minimizes XrevSE.
                scalar minVal = XrevSE[idx];
                for (label i : idx_gut)
                {
                    if (XrevSE[i] < minVal)
                    {
                        minVal = XrevSE[i];
                        idx = i;
                    }
                }
            }
        }
        else
        {
            // If no indices satisfy the criteria, choose the index that minimizes XrevSE overall.
            idx = 0;
            scalar minVal = XrevSE[0];
            for (label i = 1; i < h ; ++i)
            {
                if (XrevSE[i] < minVal)
                {
                    minVal = XrevSE[i];
                    idx = i;
                }
            }
        }

        tmin = idx*(static_cast<scalar>(k)/static_cast<scalar>(h)) + 1;
        minTimes.append(tmin);

        XrevSE = fractionalFilter(XrevSE);
        h = XrevSE.size();
    }

    if (minTimes.size() > 4)
    {
        minTimes.shrink();
        SortableList<scalar> shortTime(minTimes);
        shortTime.sort();
        scalar Q1Index = shortTime.size()/4.;
        scalar Q3Index = 3*Q1Index;
        scalar Q1median =
            (shortTime[std::floor(Q1Index)] + shortTime[std::ceil(Q1Index)])/2;
        scalar Q3median =
            (shortTime[std::floor(Q3Index)] + shortTime[std::ceil(Q3Index)])/2;

        scalar IQR = Q3median - Q1median;
        label realData = 0;
        meanTime = 0;
        forAll(shortTime, i)
        {
            if
            (
                shortTime[i] >= Q1median - 1.5*IQR
             && shortTime[i] <= Q3median + 1.5*IQR
            )
            {
                meanTime += shortTime[i];
                realData++;
            }
        }

        if (realData>0)
        {
            meanTime /= realData;
        }
        else
        {
            meanTime = sum(shortTime)/shortTime.size();
        }
    }

    // After the loop, use floor(tmin) to index into the original t vector.
    label idxT = static_cast<label>(std::floor(meanTime));

    return t[idxT];
}

Foam::List<Foam::scalar> Foam::functionObjects::runTimeControls
::confidenceIntervalCondition
::fractionalFilter
(
    const List<Foam::scalar>& X
)
{
    label N = X.size();

    label h =  static_cast<label>(std::ceil(static_cast<scalar>(N) / 2));

    List<scalar> Xnew(h, 0.0);

    if (N % 2 == 0 && N != 0)
    {
        // For even N: sum pairs of elements.
        label j = 0;
        for (label i = 0; i < h; ++i, j += 2)
        {
            Xnew[i] = X[j] + X[j + 1];
        }
    }
    else
    {
        // For odd N: use weighted sum.
        scalar w = (static_cast<scalar>(N) / h) - 1;
        label j = 0;
        for (label i = 0; i < h; ++i, j += 2)
        {
            scalar ws = static_cast<scalar>(i) / h;
            scalar we = w - ws;
            // Mimic Python negative indexing: if j-1 < 0, use the last element.
            scalar x_prev = (j - 1 < 0) ? X.last() : X[j - 1];
            Xnew[i] = ws * x_prev + X[j] + boundLimit(X, j) * we;
        }
    }
    // Multiply Xnew element-wise by (h/N)
    scalar factor = static_cast<scalar>(h) / N;
    for (label i = 0; i < h; ++i)
    {
        Xnew[i] *= factor;
    }
    return Xnew;
}

Foam::scalar Foam::functionObjects::runTimeControls::confidenceIntervalCondition
::boundLimit
(
    const List<scalar>& X, label i
)
{
    if (i + 1 < X.size())
    {
        return X[i + 1];
    }
    else
    {
        return X.last();
    }
}

Foam::List<Foam::scalar> Foam::functionObjects::runTimeControls
::confidenceIntervalCondition::reverseStandardError
(
    const List<scalar>& X, label h
)
{
    List<scalar> XrevSE(h, 0.0);

    scalar mV = 0;
    scalar oldmV = mV;
    scalar M2 = 0;

    for (label j = 0; j < h; ++j)
    {
        scalar beta = 1./(j+1);
        if (j==0)
        {
            mV = beta*X[h-1-j];
            XrevSE[h-1-j] = 0;
        }
        else
        {
            mV = (1-beta)*oldmV + beta*X[h-1-j];
            M2 = (1-beta)*M2 + beta*(X[h-1-j]-oldmV)*(X[h-1-j]-mV);
            XrevSE[h-1-j] = std::sqrt(M2/(j+1));
        };

        oldmV = mV;
    }

    auto cutOff =  static_cast<label>(std::ceil(0.1*static_cast<scalar>(h)));
    for (label j = h-cutOff; j < h; ++j)
    {
        XrevSE[j] = XrevSE[h-cutOff];
    }

    return XrevSE;
}

Foam::scalar Foam::functionObjects::runTimeControls
::confidenceIntervalCondition::meanData
(
    const List<Foam::scalar>& X, Foam::label start, Foam::label end
)
{
    scalar sum =0;
    for (int i = start; i<end; ++i)
    {
        sum += X[i];
    }
    return sum / (end - start);
}

Foam::scalar Foam::functionObjects::runTimeControls
::confidenceIntervalCondition::stErr
(
    const List<Foam::scalar>& X, Foam::label start, Foam::label end, scalar mean
)
{
    scalar sqSum = 0.0;


    for (label i = start; i < end; ++i)
    {
        scalar diff = X[i] - mean;
        sqSum += sqr(diff);
    }

    scalar n = end - start;
    scalar stddev = std::sqrt(sqSum / n);
    return stddev / std::sqrt(n);
}

Foam::List<Foam::scalar> Foam::functionObjects::runTimeControls
::confidenceIntervalCondition::autoCorrelationFactor
(
    const List<Foam::scalar>& X
)
{
    label N = X.size();
    scalar meanX = sum(X)/N;
    scalarField Xstd(X-meanX);
    DynamicList<scalar> autocorrelation(0);

    scalar dataNorm = 0;

    forAll(Xstd, i)
    {
        dataNorm += sqr(Xstd[i]);
    }

    autocorrelation.append(1);
    for (int i=1; i<N; i++)
    {
        scalar sumProduct = 0;
        label end = N - i;
        for (int j=0; j<end; j++)
        {
            sumProduct += Xstd[j] * Xstd[i+j];
        }

        if (dataNorm>ROOTVSMALL)
        {
            sumProduct /= dataNorm;
        }

        autocorrelation.append(sumProduct);
        if (autocorrelation[i]<0)
        {
            break;
        }
    }

    if (debug)
    {
        word name("autoCorrelation");
        fileName fileN(name);
        OFstream os(fileN);
        forAll(autocorrelation, i)
        {
            os << autocorrelation[i] << nl;
        }
    }

    DynamicList<scalar> regressedAutoCorr(0);

    forAll(autocorrelation, i)
    {
        if (autocorrelation[i]>0)
        {
            regressedAutoCorr.append(autocorrelation[i]);
        }
        else
        {
            break;
        }
    }

    return regressedAutoCorr.shrink();
}

Foam::scalar Foam::functionObjects::runTimeControls
::confidenceIntervalCondition::autoCorrelationTime
(
    const List<Foam::scalar>& autoCorrFactor
)
{
    scalarField xAc(autoCorrFactor.size());
    for (label i = 0; i < autoCorrFactor.size(); i++)
    {
        xAc[i] = scalar(i);
    }

    // Fit an exponential decay to autoCorrFactor
    // using linear regression on ln(acf)
    scalar tau = simpleLinearRegression(xAc, autoCorrFactor);

    // Compute autocorrelation time as the half-life:
    //    autocorrelation time = ceil(tau * ln2)
    scalar autoTime = std::ceil(tau * std::log(2.0));

    return autoTime;
}

Foam::scalar Foam::functionObjects::runTimeControls
::confidenceIntervalCondition::simpleLinearRegression
(
    const List<Foam::scalar>& X, const List<Foam::scalar>& Y
)
{
    if (X.size()==1)
    {
        return 1;
    }
    scalar sumX = 0.0;
    scalar sumY = 0.0;
    scalar sumXY = 0.0;
    scalar sumX2 = 0.0;
    const label n = X.size();

    for (label i = 0; i < n; i++)
    {
        // Assuming acf[i] > 0; otherwise, log() is undefined.
        const scalar xi = X[i];
        const scalar yi = std::log(Y[i]);

        sumX += xi;
        sumY += yi;
        sumXY += xi * yi;
        sumX2 += xi * xi;
    }

    const scalar slope = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);

    return -1.0 / slope;
}

void Foam::functionObjects::runTimeControls
::confidenceIntervalCondition::computeEquilibratedUncertainty
(
    const List<Foam::scalar>& X,
    const Foam::scalar& aCorrTime,
    scalar& equilibratedAverage,
    scalar& equilibratedUncertainty,
    label& uncorrelatedBatches
)
{
    scalarField uncorrBatchesValues(0);
    // Ensure acTime is an integer by applying ceiling.
    const auto acTimeInt = static_cast<label>(std::ceil(aCorrTime));

    // Divide data into uncorrelated chunks if acTimeInt > 1
    if (acTimeInt > 1)
    {
        // Determine the number of complete batches.
        const label nBatches = X.size() / acTimeInt;

        // Allocate uncorrBatches with nBatches elements.
        uncorrBatchesValues.setSize(nBatches);

        // Compute the average for each batch.
        for (label i = 0; i < nBatches; i++)
        {
            scalar sum = 0.0;
            for (label j = 0; j < acTimeInt; j++)
            {
                sum += X[i * acTimeInt + j];
            }
            uncorrBatchesValues[i] = sum / acTimeInt;
        }
    }
    else
    {
        // When acTimeInt is 1, each data point is treated as an uncorrelated batch.
        uncorrBatchesValues = X;
    }

    equilibratedAverage = sum(X)/X.size();

    uncorrelatedBatches = uncorrBatchesValues.size();
    scalar meanBatch = sum(uncorrBatchesValues)/uncorrelatedBatches;

    scalar sqSum = 0.0;
    for (label i = 0; i < uncorrelatedBatches; i++)
    {
        const scalar diff = uncorrBatchesValues[i] - meanBatch;
        sqSum += diff * diff;
    }
    const scalar stdDev = sqrt(sqSum/uncorrelatedBatches);

    // Standard error of the mean.
    equilibratedUncertainty = stdDev / sqrt(static_cast<scalar>(uncorrelatedBatches));
}

Foam::scalar Foam::functionObjects::runTimeControls::confidenceIntervalCondition
::getRegressedDataSlope
(
    const List<scalar>& T,
    const List<scalar>& X
)
{
    if (X.size()==1)
    {
        return 0;
    }

    scalar sumX = 0.0;
    scalar sumY = 0.0;
    scalar sumXY = 0.0;
    scalar sumX2 = 0.0;
    const label n = T.size();

    for (label i = 0; i < n; i++)
    {
        // Assuming acf[i] > 0; otherwise, log() is undefined.
        const scalar xi = T[i];
        const scalar yi = X[i];

        sumX += xi;
        sumY += yi;
        sumXY += xi * yi;
        sumX2 += xi * xi;
    }

    scalar slope = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);

    return slope;
}


bool Foam::functionObjects::runTimeControls::confidenceIntervalCondition
::checkConvergence
(
    scalar tolerance,
    scalar confidenceInterval,
    scalar slope,
    scalar size
) const
{
    if
    (
        mag(slope) < slopeTolerance_
     && confidenceInterval < tolerance
     && size > minSamples_
    )
    {
        return true;
    }

    return  false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls
::confidenceIntervalCondition::confidenceIntervalCondition
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    stateFunctionObject& state
)
:
    runTimeCondition(name, obr, dict, state),
    writeFile(obr_, name, typeName, dict),
    dict_(dict),
    functionObjectName_(dict.lookupOrDefault<word>("functionObject","none")),
    entryName_(dict.lookup("field")),
    confidence_(dict.lookup<scalar>("confidence")),
    interval_(dict.lookup<scalar>("targetInterval")),
    intervalType_
    (
        intvTypeNames_.read
        (
            IStringStream
            (
                dict.lookupOrDefault<word>
                (
                    "intervalType",
                    word("relative")
                )
            )()
        )
    ),
    minSamples_(dict.lookupOrDefault<label>("minSamples", 500)),
    tDistFile_
    (
        dict.lookupOrDefault<fileName>
        (
            "distributionFile",
            "dictData/tDistribution.dat"
        )
    ),
    filterRadius_(dict.lookupOrDefault<label>("filterRadius", 0)),
    slopeTolerance_(dict.lookupOrDefault<scalar>("slopeTolerance", 5e-3)),
    restart_(dict.lookupOrDefault<Switch>("restart", true)),
    smoothValueForTransientDetection_
    (
        dict.lookupOrDefault<Switch>("smoothForTransientDetection", false)
    ),
    smoothValueForConvergence_
    (
        dict.lookupOrDefault<Switch>("smoothForConvergence", true)
    ),
    minTransientTime_(dict.lookupOrDefault<scalar>("minTransientTime", 0))

{
    IFstream distFile(findEtcFile(tDistFile_));

    dictionary tDict(distFile);
    tDistribution_ = Function2<scalar>::New("distribution", tDict);

    word previousPath = fileName(conditionDict_.lookupOrDefault<fileName>
        ("previousSamplingData", word::null));

    if (previousPath != word::null && restart_)
    {
        IFstream is(previousPath);
        string line;
        // Skip Header
        is.getLine(line);
        while (is.getLine(line))
        {
            std::istringstream iss(line);
            scalar nSamples, time, transientTime, value, mean;
            iss >> nSamples >> time >> transientTime >> value >> mean;

            if (time <= obr_.time().value())
            {
                T_.append(time);
                X_.append(value);
            }
        }
    }

    //set up file IO
    if (Pstream::master() || !Pstream::parRun())
    {
        //add headers
        writeCommented(file(), "nSamples");
        writeDelimited(file(), "Time");
        writeDelimited(file(), "Transient_Time");
        writeDelimited(file(), "Value");
        writeDelimited(file(), "Mean");
        writeDelimited(file(), "Upper");
        writeDelimited(file(), "Lower");
        writeDelimited(file(), "Filtered");
        writeDelimited(file(), "Slope");

        file() <<endl;
    }

}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls
::confidenceIntervalCondition::~confidenceIntervalCondition()
= default;


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

bool Foam::functionObjects::runTimeControls
::confidenceIntervalCondition::apply(bool postProcess)
{
    bool satisfied = true;

    if (!active_)
    {
        return satisfied;
    }

    if (postProcess)
    {
        satisfied = false;

        postProcessMode(satisfied);

        return satisfied;
    }

    scalar t = obr_.time().value();

    satisfied = false;

    Log << "    " << type() << ": ";

    auto X
        = state_.getObjectResult<scalar>(functionObjectName_, entryName_);

    X_.append(X);
    T_.append(t);

    if (X_.size() > filterRadius_)
    {
        calc(satisfied);
    }

    Log << endl;

    return satisfied;
}

void Foam::functionObjects::runTimeControls::confidenceIntervalCondition
::postProcessMode(bool& satisfied)
{
    word pathName = dict_.lookupType<word>("dataFile");
    IFstream is(pathName);
    string line;

    // Skip Header
    is.getLine(line);
    std::istringstream firstLine(line);
    DynamicList<word> listOfFields;
    word fieldName;

    while (firstLine >> fieldName)
    {
        listOfFields.append(fieldName);
    }
    listOfFields.shrink();
    bool foundField = false;
    label fieldIndex = -1;
    forAll(listOfFields, i)
    {
        if (listOfFields[i] == entryName_)
        {
            foundField = true;
            fieldIndex = i;
        }
    }

    if (!foundField)
    {
        FatalErrorInFunction
        << "Field "<<entryName_
        <<" is not in the file you specified "
        << exit(FatalError);
    }

    while (is.getLine(line))
    {
        std::istringstream iss(line);
        scalar time, X;
        label counter = 0;
        iss >> time;

        while (iss >> X)
        {
            if (counter == fieldIndex-2)
            {
                break;
            }
            counter++;
        }

        X_.append(X);
        T_.append(time);

        if (X_.size() > filterRadius_)
        {
            calc(satisfied);
        }

        if (satisfied)
        {
            return;
        }
    }
}


void Foam::functionObjects::runTimeControls
::confidenceIntervalCondition::write()
{

    fileName lastDataSeries
    (
        "postProcessing/"
      + name_ + "/"
      + obr_.time().timeName(obr_.time().startTime().value()) + "/"
      + typeName + ".dat"
    );

    conditionDict_.add("previousSamplingData", lastDataSeries);

}


// ************************************************************************* //
