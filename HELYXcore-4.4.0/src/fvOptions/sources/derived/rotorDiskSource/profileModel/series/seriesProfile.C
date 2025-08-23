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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "seriesProfile.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/IOstreams/Fstreams/IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(seriesProfile, 0);
    addToRunTimeSelectionTable(profileModel, seriesProfile, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::seriesProfile::evaluateDrag
(
    const scalar& xIn,
    const List<scalar>& values
) const
{
    scalar result = 0.0;

    forAll(values, i)
    {
        result += values[i]*cos(i*xIn);
    }

    return result;
}


Foam::scalar Foam::seriesProfile::evaluateLift
(
    const scalar& xIn,
    const List<scalar>& values
) const
{
    scalar result = 0.0;

    forAll(values, i)
    {
        // note: first contribution always zero since sin(0) = 0, but
        // keep zero base to be consitent with drag coeffs
        result += values[i]*sin(i*xIn);
    }

    return result;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::seriesProfile::seriesProfile
(
    const dictionary& dict,
    const word& modelName
)
:
    profileModel(dict, modelName),
    CdCoeffs_(),
    ClCoeffs_()
{
    if (readFromFile())
    {
        IFstream is(fName_);
        is  >> CdCoeffs_ >> ClCoeffs_;
    }
    else
    {
        CdCoeffs_ = dict.lookup<List<scalar>>("CdCoeffs");
        ClCoeffs_ = dict.lookup<List<scalar>>("ClCoeffs");
    }


    if (CdCoeffs_.empty())
    {
        FatalErrorInFunction
            << "CdCoeffs must be specified" << exit(FatalError);
    }
    if (ClCoeffs_.empty())
    {
        FatalErrorInFunction
            << "ClCoeffs must be specified" << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::seriesProfile::Cdl(const scalar alpha, scalar& Cd, scalar& Cl) const
{
    Cd = evaluateDrag(alpha, CdCoeffs_);
    Cl = evaluateLift(alpha, ClCoeffs_);
}


// ************************************************************************* //
