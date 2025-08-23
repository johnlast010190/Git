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
    (c) 2013-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "turbulentBreakUp.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace IATEsources
{
    defineTypeNameAndDebug(turbulentBreakUp, 0);
    addToRunTimeSelectionTable(IATEsource, turbulentBreakUp, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::IATEsources::turbulentBreakUp::
turbulentBreakUp
(
    const IATE& iate,
    const dictionary& dict
)
:
    IATEsource(iate),
    Cti_("Cti", dimless, dict),
    WeCr_("WeCr", dimless, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::diameterModels::IATEsources::turbulentBreakUp::R() const
{
    tmp<volScalarField> tR
    (
        volScalarField::New
        (
            "R",
            iate_.phase().U().mesh(),
            dimensionedScalar(dimless/dimTime, 0)
        )
    );

    volScalarField R = tR();

    scalar Cti = Cti_.value();
    scalar WeCr = WeCr_.value();
    volScalarField Ut(this->Ut());
    volScalarField We(this->We());

    tmp<volScalarField> td(iate_.d());
    const volScalarField& d = td();

    forAll(R, celli)
    {
        if (We[celli] > WeCr)
        {
            R[celli] =
                (1.0/3.0)*Cti/d[celli]*Ut[celli]*sqrt(1 - WeCr/We[celli])
               *exp(-WeCr/We[celli]);
        }
    }

    return tR;
}


// ************************************************************************* //
