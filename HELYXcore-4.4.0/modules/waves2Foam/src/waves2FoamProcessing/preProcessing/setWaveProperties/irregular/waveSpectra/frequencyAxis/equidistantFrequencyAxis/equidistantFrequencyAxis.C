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

#include "equidistantFrequencyAxis.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(equidistantFrequencyAxis, 0);
addToRunTimeSelectionTable(frequencyAxis, equidistantFrequencyAxis, frequencyAxis);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


equidistantFrequencyAxis::equidistantFrequencyAxis
(
    const Time& rT,
    dictionary& dict
)
:
    frequencyAxis(rT, dict)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalarField equidistantFrequencyAxis::freqAxis
(
    const scalarField&,
    const scalarField&,
    const label& N
) const
{
    scalarField freq(N + 1, 0.0);

    scalar df = (fu_ - fl_)/static_cast<scalar>(N);

    for (int i = 0; i < N + 1; i++)
    {
        freq[i] = static_cast<scalar>(i)*df + fl_;
    }

    return freq;
}


scalarField equidistantFrequencyAxis::freqAxis(const label& N) const
{
    scalarField dummy(0);

    return freqAxis(dummy, dummy, N);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
