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
    (c) held by original author

\*---------------------------------------------------------------------------*/

#include "frequencyAxis.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(frequencyAxis, 0);
defineRunTimeSelectionTable(frequencyAxis, frequencyAxis);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


frequencyAxis::frequencyAxis
(
    const Time& rT,
    dictionary& dict
)
:
    rT_(rT),
    dict_(dict)
{
    if (dict_.subDict("frequencyAxis").found("lowerFrequencyCutoff"))
    {
        fl_ = readScalar
            (
                dict_.subDict("frequencyAxis").lookup("lowerFrequencyCutoff")
            );
    }
    else
    {
        scalar Tp = dict_.lookup<scalar>("Tp");
        fl_ = 0.3/Tp;
    }

    if (dict_.subDict("frequencyAxis").found("upperFrequencyCutoff"))
    {
        fu_ = readScalar
            (
                dict_.subDict("frequencyAxis").lookup("upperFrequencyCutoff")
            );
    }
    else
    {
        scalar Tp = dict_.lookup<scalar>("Tp");
        fu_ = 3.0/Tp;
    }
}


frequencyAxis::~frequencyAxis()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


autoPtr<frequencyAxis> frequencyAxis::New
(
    const Time& rT,
    dictionary& dict
)
{
    word discretisation = dict.subDict("frequencyAxis").lookup("discretisation");

    const auto ctor =
        ctorTableLookup
        (
            "discretisation method",
            frequencyAxisConstructorTable_(),
            discretisation
        );
    return autoPtr<frequencyAxis>(ctor(rT, dict));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
