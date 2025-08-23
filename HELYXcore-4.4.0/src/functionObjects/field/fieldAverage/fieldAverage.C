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
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fieldAverage/fieldAverage.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldAverage, 0);
    addToRunTimeSelectionTable(functionObject, fieldAverage, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldAverage::initialise()
{
    // Initialise any unset times
    forAll(totalTime_, fieldi)
    {
        if (totalTime_[fieldi] < 0)
        {
            totalTime_[fieldi] = obr().time().deltaTValue();
        }
    }

    // Initialise mean fields
    forAll(faItems_, fieldi)
    {
        initialiseMeanField<scalar>(fieldi);
        initialiseMeanField<vector>(fieldi);
        initialiseMeanField<sphericalTensor>(fieldi);
        initialiseMeanField<symmTensor>(fieldi);
        initialiseMeanField<tensor>(fieldi);
    }

    // Initialise prime-squared mean fields
    forAll(faItems_, fieldi)
    {
        initialisePrime2MeanField<scalar, scalar>(fieldi);

        if (!Welford_)
        {
            initialisePrime2MeanField<vector, symmTensor>(fieldi);
        }
        else
        {
            if (!symmetrize_)
            {
                initialisePrime2MeanFieldWel<vector, tensor>(fieldi);
            }
            else
            {
                initialisePrime2MeanFieldWelSymm<vector, symmTensor>(fieldi);
            }
        }
    }

    // Add min fields
    forAll(faItems_, fieldi)
    {
        addMinField<scalar>(fieldi);
        addMinField<vector>(fieldi);
        addMinField<sphericalTensor>(fieldi);
        addMinField<symmTensor>(fieldi);
        addMinField<tensor>(fieldi);
    }

    // Add max fields
    forAll(faItems_, fieldi)
    {
        addMaxField<scalar>(fieldi);
        addMaxField<vector>(fieldi);
        addMaxField<sphericalTensor>(fieldi);
        addMaxField<symmTensor>(fieldi);
        addMaxField<tensor>(fieldi);
    }
}


void Foam::functionObjects::fieldAverage::restart()
{
    Log << "    Restarting averaging at time " << obr().time().timeName() << nl;

    // Clear the times
    totalIter_ = 1;
    totalTime_ = -1;

    // Clear mean fields
    forAll(faItems_, fieldi)
    {
        if (faItems_[fieldi].mean())
        {
            if (obr().found(faItems_[fieldi].meanFieldName()))
            {
                obr().checkOut(*obr()[faItems_[fieldi].meanFieldName()]);
            }
        }
        if (faItems_[fieldi].prime2Mean())
        {
            if (obr().found(faItems_[fieldi].prime2MeanFieldName()))
            {
                obr().checkOut(*obr()[faItems_[fieldi].prime2MeanFieldName()]);
            }
        }
        if (faItems_[fieldi].minFlag())
        {
            if (obr().found(faItems_[fieldi].minFieldName()))
            {
                obr().checkOut(*obr()[faItems_[fieldi].minFieldName()]);
            }
        }
        if (faItems_[fieldi].maxFlag())
        {
            if (obr().found(faItems_[fieldi].maxFieldName()))
            {
                obr().checkOut(*obr()[faItems_[fieldi].maxFieldName()]);
            }
        }
    }

    // Re-create any mean fields
    initialise();

    // Ensure first averaging works unconditionally
    prevTimeIndex_ = -1;
}


void Foam::functionObjects::fieldAverage::calcAverages()
{
    Log << type() << " " << name() << ":" << nl;

    const label currentTimeIndex = obr().time().timeIndex();
    const scalar currentTime = obr().time().value();

    if (prevTimeIndex_ == currentTimeIndex)
    {
        return;
    }

    prevTimeIndex_ = currentTimeIndex;

    if (periodicRestart_ && currentTime > restartPeriod_*periodIndex_)
    {
        restart();
        periodIndex_++;
    }
    else
    {
        initialise();
    }

    Log << "    Calculating averages" << nl;

    if (!Welford_)
    {
        addMeanSqrToPrime2Mean<scalar, scalar>();
        addMeanSqrToPrime2Mean<vector, symmTensor>();

        calculateMeanFields<scalar>();
        calculateMeanFields<vector>();

        calculatePrime2MeanFields<scalar, scalar>();
        calculatePrime2MeanFields<vector, symmTensor>();
    }
    else
    {
        calculateMeanPrime2MeanFieldsWel<scalar, scalar>();
        if (!symmetrize_)
        {
            calculateMeanPrime2MeanFieldsWel<vector, tensor>();
        }
        else
        {
            calculateMeanPrime2MeanFieldsWelSymm<vector, symmTensor>();
        }
    }
    calculateMeanFields<sphericalTensor>();
    calculateMeanFields<symmTensor>();
    calculateMeanFields<tensor>();

    calculateMinFields<scalar>();
    calculateMinFields<vector>();
    calculateMinFields<sphericalTensor>();
    calculateMinFields<symmTensor>();
    calculateMinFields<tensor>();

    calculateMaxFields<scalar>();
    calculateMaxFields<vector>();
    calculateMaxFields<sphericalTensor>();
    calculateMaxFields<symmTensor>();
    calculateMaxFields<tensor>();

    forAll(faItems_, fieldi)
    {
        totalIter_[fieldi]++;
        totalTime_[fieldi] += obr().time().deltaTValue();
    }

    Log << endl;
}


void Foam::functionObjects::fieldAverage::writeAverages()
{
    Log << "    Writing average fields" << nl << endl;

    writeFields<scalar>();
    writeFields<vector>();
    writeFields<sphericalTensor>();
    writeFields<symmTensor>();
    writeFields<tensor>();

    forAll(faItems_, fieldi)
    {
        const word& fieldName = faItems_[fieldi].fieldName();

        dictionary propsDict;
        propsDict.add("totalIter", totalIter_[fieldi]);
        propsDict.add("totalTime", totalTime_[fieldi]);
        setProperty(fieldName, propsDict);
    }
}


void Foam::functionObjects::fieldAverage::read
(
    const dictionary& dict,
    const bool construct
)
{
    const scalar currentTime = obr().time().value();

    dict.readIfPresent("restartOnRestart", restartOnRestart_);
    dict.readIfPresent("restartOnOutput",  restartOnOutput_);
    dict.readIfPresent("periodicRestart",  periodicRestart_);

    if (periodicRestart_)
    {
        scalar userRestartPeriod = dict.lookup<scalar>("restartPeriod");
        restartPeriod_ = obr().time().userTimeToTime(userRestartPeriod);

        if (restartPeriod_ > 0)
        {
            // Determine the appropriate interval for the next restart
            periodIndex_ = 1;
            while (currentTime > restartPeriod_*periodIndex_)
            {
                ++periodIndex_;
            }

            Info<< "    Restart period " << userRestartPeriod
                << " - next restart at " << (userRestartPeriod*periodIndex_)
                << nl << endl;
        }
        else
        {
            periodicRestart_ = false;

            Info<< "    Restart period " << userRestartPeriod
                << " - ignored" << nl << endl;
        }
    }

    scalar userRestartTime = 0;
    if (dict.readIfPresent("restartTime", userRestartTime))
    {
        restartTime_ = obr().time().userTimeToTime(userRestartTime);

        if (currentTime > restartTime_)
        {
            // The restart time is already in the past - ignore.
            restartTime_ = GREAT;
        }
        else
        {
            Info<< "    Restart scheduled at time " << userRestartTime
                << nl << endl;
        }
    }

    if (dict.found("averagingStartTime"))
    {
        Warning << "Use of 'averagingStartTime' keyword in fieldAverage "
                << "functionObject is deprecated." << nl << "Use the general "
                << "functionObject keyword 'timeStart' instead." << endl;
    }

    dict.readIfPresent("Welford", Welford_);
    dict.readIfPresent("symmetrize", symmetrize_);

    if (construct)
    {
        // First read of a run

        faItems_.clear();
        totalIter_.clear();
        totalTime_.clear();

        // Read fields
        faItems_ = dict.lookup<List<fieldAverageItem>>("fields");

        totalIter_.setSize(faItems_.size(), 1);
        totalTime_.setSize(faItems_.size(), -1);

        bool first = true;
        forAll(faItems_, fieldi)
        {
            const word& fieldName = faItems_[fieldi].fieldName();
            if (!foundProperty(fieldName))
            {
                if (first)
                {
                    Log << "    Starting averaging for fields:" << nl;
                    first = false;
                }
                Log << "        " << fieldName << nl;
            }
        }

        first = true;
        forAll(faItems_, fieldi)
        {
            const word& fieldName = faItems_[fieldi].fieldName();
            if (foundProperty(fieldName))
            {
                if (first)
                {
                    Log << "    Restarting averaging for fields:" << nl;
                    first = false;
                }

                dictionary fieldDict;
                getDict(fieldName, fieldDict);

                totalIter_[fieldi] = fieldDict.lookup<label>("totalIter");
                totalTime_[fieldi] = fieldDict.lookup<scalar>("totalTime");

                scalar userTotalTime =
                    obr().time().timeToUserTime(totalTime_[fieldi]);

                Info<< "        " << fieldName
                    << " iters = " << totalIter_[fieldi]
                    << " time = " << userTotalTime << nl;
            }
        }
    }
    else
    {
        // Re-read during a run. Read the items and copy the per-field total
        // iter/times from the old items to the new.

        List<fieldAverageItem> faItems0;
        List<label> totalIter0;
        List<scalar> totalTime0;
        faItems0.transfer(faItems_);
        totalIter0.transfer(totalIter_);
        totalTime0.transfer(totalTime_);

        faItems_ = dict.lookup<List<fieldAverageItem>>("fields");

        totalIter_.setSize(faItems_.size(), 1);
        totalTime_.setSize(faItems_.size(), -1);

        // Map from field to old-field index
        labelList fieldiFieldi0s(faItems_.size(), -1);
        forAll(faItems_, fieldi)
        {
            const word& fieldName = faItems_[fieldi].fieldName();
            forAll(faItems0, fieldi0)
            {
                if (faItems0[fieldi0].fieldName() == fieldName)
                {
                    fieldiFieldi0s[fieldi] = fieldi0;
                    break;
                }
            }
        }

        bool first = true;
        forAll(faItems_, fieldi)
        {
            if (fieldiFieldi0s[fieldi] == -1)
            {
                if (first)
                {
                    Log << "    Starting averaging for fields:" << nl;
                    first = true;
                }
                Log << "        " << faItems_[fieldi].fieldName() << nl;
            }
        }

        first = true;
        forAll(faItems_, fieldi)
        {
            if (fieldiFieldi0s[fieldi] != -1)
            {
                if (first)
                {
                    Log << "    Continuing averaging for fields:" << nl;
                    first = true;
                }
                totalIter_[fieldi] = totalIter0[fieldiFieldi0s[fieldi]];
                totalTime_[fieldi] = totalTime0[fieldiFieldi0s[fieldi]];
                Log << "        " << faItems_[fieldi].fieldName()
                    << " iters = " << totalIter_[fieldi]
                    << " time = " << totalTime_[fieldi] << nl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldAverage::fieldAverage
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    prevTimeIndex_(-1),
    restartOnRestart_(false),
    restartOnOutput_(false),
    periodicRestart_(false),
    restartPeriod_(GREAT),
    restartTime_(GREAT),
    Welford_(false),
    symmetrize_(true),
    faItems_(),
    totalIter_(),
    totalTime_(),
    periodIndex_(1)
{
    fvMeshFunctionObject::read(dict);

    Log << type() << " " << name << ":" << nl;

    read(dict, true);

    // Read any available mean fields
    forAll(faItems_, fieldi)
    {
        readMeanField<scalar>(fieldi);
        readMeanField<vector>(fieldi);
        readMeanField<sphericalTensor>(fieldi);
        readMeanField<symmTensor>(fieldi);
        readMeanField<tensor>(fieldi);
    }

    // Read any available prime-squared mean fields
    forAll(faItems_, fieldi)
    {
        readPrime2MeanField<scalar, scalar>(fieldi);
        readPrime2MeanField<vector, symmTensor>(fieldi);
    }

    Log << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldAverage::~fieldAverage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldAverage::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Log << type() << " " << name() << ":" << nl;

    read(dict, false);

    Log << endl;

    return true;
}


bool Foam::functionObjects::fieldAverage::execute()
{
    calcAverages();

    return true;
}


bool Foam::functionObjects::fieldAverage::write()
{
    writeAverages();

    if (restartOnOutput_)
    {
        restart();
    }

    return true;
}


// ************************************************************************* //
