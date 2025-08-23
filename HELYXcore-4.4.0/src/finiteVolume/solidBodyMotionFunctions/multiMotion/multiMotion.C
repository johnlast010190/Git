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
    (c) 2019-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "solidBodyMotionFunctions/multiMotion/multiMotion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(multiMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        multiMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::multiMotion::multiMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime,
    const word& frameName
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime, frameName)
{
    DeprecationWarningInFunction
    (
        "multiMotion",
        "motion type",
        40200,
        "Please replace it with nested dynamic frames."
    );

    multiMotion::read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::multiMotion::~multiMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::multiMotion::transformation
(
    const scalar t1,
    const scalar t2
) const
{
    const septernion TR = multiMotionTransform(0.0, t2);
    if (isIncrementalMotion())
    {
        // check all motions
        for (label i = 0; i < SBMFs_.size(); i++)
        {
            if (SBMFs_[i].isIncrementalMotion())
            {
                FatalErrorInFunction
                    << "Inconsistent multiMotion setup. Associated motion"
                    << nl << "functions must all be absolute."
                    << exit(FatalError);
            }
        }

        return TR/multiMotionTransform(0.0, t2 - time_.deltaT().value());
    }

    DebugInFunction
        << "Transformation from time " << t1
        << " to " << t2 << " is " << TR << endl;

    return TR;
}


Foam::septernion
Foam::solidBodyMotionFunctions::multiMotion::multiMotionTransform
(
    const scalar t1,
    const scalar t2
) const
{
    septernion TR = SBMFs_[0].transformation(t1, t2);

    for (label i = 1; i < SBMFs_.size(); i++)
    {
        TR *= SBMFs_[i].transformation(t1, t2);
    }

    return TR;
}


bool Foam::solidBodyMotionFunctions::multiMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    label i = 0;
    SBMFs_.setSize(SBMFCoeffs_.size());

    forAllConstIter(IDLList<entry>, SBMFCoeffs_, iter)
    {
        if (iter().isDict())
        {
            SBMFs_.set(i, solidBodyMotionFunction::New(iter().dict(), time_));

            Info<< "Constructed SBMF " << i << " : " << iter().keyword()
                << " of type " << SBMFs_[i].type() << endl;
            i++;
        }
    }
    SBMFs_.setSize(i);

    return true;
}


// ************************************************************************* //
