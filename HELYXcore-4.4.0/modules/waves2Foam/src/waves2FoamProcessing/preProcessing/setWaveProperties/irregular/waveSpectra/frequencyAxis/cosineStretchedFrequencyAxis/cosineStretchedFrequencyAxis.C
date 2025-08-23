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

#include "cosineStretchedFrequencyAxis.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(cosineStretchedFrequencyAxis, 0);
addToRunTimeSelectionTable(frequencyAxis, cosineStretchedFrequencyAxis, frequencyAxis);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


cosineStretchedFrequencyAxis::cosineStretchedFrequencyAxis
(
    const Time& rT,
    dictionary& dict
)
:
    frequencyAxis(rT, dict)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalarField cosineStretchedFrequencyAxis::freqAxis
(
    const scalarField&,
    const scalarField&,
    const label& N
) const
{
    scalarField freq(N + 1, 0.0);

    scalar fp = 1.0/dict_.lookup<scalar>("Tp");

    // Calculate the number of upper and lower frequencies
    label Nlow(ceil((fp - fl_)/(fu_ - fp)*(N + 1)));
    label Nhigh(N + 1 - Nlow);

    for (int i = 0; i < Nlow; i++)
    {
        freq[i] = (fp - fl_)*Foam::sin(2*M_PI/(4.0*Nlow)*i) + fl_;
    }

    for (int i = 0; i <= Nhigh; i++)
    {
        freq[Nlow - 1 + i] =
            (fu_ - fp)*(- Foam::cos(2*M_PI/(4*Nhigh)*i) + 1) + fp;
    }

    return freq;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
