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

#include "waveSpectra.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(waveSpectra, 0);
defineRunTimeSelectionTable(waveSpectra, waveSpectra);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


waveSpectra::waveSpectra
(
    const Time& rT,
    dictionary& dict,
    scalarField& amp,
    scalarField& freq,
    scalarField& phi,
    vectorField& k
)
:
    rT_(rT),
    dict_(dict),
    amp_(amp),
    freq_(freq),
    phi_(phi),
    k_(k),

    G_
    (
        Foam::mag
        (
            uniformDimensionedVectorField
            (
                rT_.db().lookupObject<uniformDimensionedVectorField>("g")
            ).value()
        )
    ),

    PI_( M_PI ),

    phases_(Foam::phases::New(rT_, dict_))
{
}


waveSpectra::~waveSpectra()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void waveSpectra::writeSpectrum
(
    Ostream& os,
    const scalarField& freq,
    const scalarField& S
) const
{
    if (dict_.subDict("frequencyAxis").lookupOrDefault<Switch>("writeSpectrum",false))
    {
        S.writeEntry("spectramValues", os);
        os << nl;

        freq.writeEntry("fspectrum", os);
        os << nl;
    }
}


autoPtr<waveSpectra> waveSpectra::New
(
    const Time& rT,
    dictionary& dict,
    scalarField& amp,
    scalarField& freq,
    scalarField& phi,
    vectorField& k
)
{
    word spectrumName = dict.lookup<word>("spectrum");

    const auto ctor =
        ctorTableLookup
        (
            "wave spectrum",
            waveSpectraConstructorTable_(),
            spectrumName
        );
    return autoPtr<waveSpectra>(ctor(rT, dict, amp, freq, phi, k));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
