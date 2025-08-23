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
    (c) 2016 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "windowModels/Hanning/Hanning.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/constants/mathematical/mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace windowModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Hanning, 0);
addToRunTimeSelectionTable(windowModel, Hanning, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Hanning::Hanning(const dictionary& dict, const label nSamples)
:
    windowModel(dict, nSamples),
    symmetric_(dict.lookup<bool>("symmetric")),
    extended_(dict.lookup<bool>("extended")),
    alpha_(dict.lookupOrDefault("alpha", 0.5)) // Hamming = 0.54
{
    // Extend range if required
    label offset = extended_ ? 1 : 0;
    scalar m = nSamples - 1 + 2*offset;

    scalarField t(nSamples);
    forAll(t, i)
    {
        t[i] = i + offset;
    }

    scalarField& wf = *this;
    wf = alpha_ - (1 - alpha_)*cos(constant::mathematical::twoPi*t/m);

    // Reset second half of window if symmetric
    if (symmetric_)
    {
        label midPointI = 0;
        if (nSamples % 2 == 0)
        {
            midPointI = nSamples/2;
        }
        else
        {
            midPointI = (nSamples + 1)/2;
        }

        for (label i = 0; i < midPointI; i++)
        {
            wf[nSamples - i - 1] = wf[i];
        }
    }

    scalar sumSqr = sum(sqr(wf));

    // Normalisation
    wf *= sqrt(nSamples/sumSqr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Hanning::~Hanning()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Hanning::symmetric() const
{
    return symmetric_;
}


bool Hanning::extended() const
{
    return extended_;
}


Foam::scalar Hanning::alpha() const
{
    return alpha_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace windowModels
} // End namespace Foam

// ************************************************************************* //
