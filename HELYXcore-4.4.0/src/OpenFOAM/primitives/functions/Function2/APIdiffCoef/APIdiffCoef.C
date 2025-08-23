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
    (c) 2011-2020 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function2/APIdiffCoef/APIdiffCoef.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function2s
{
    makeScalarFunction2(APIdiffCoef)
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function2s::APIdiffCoef::APIdiffCoef
(
    const word& name,
    const scalar a,
    const scalar b,
    const scalar wf,
    const scalar wa
)
:
    FieldFunction2<scalar, APIdiffCoef>(name),
    a_(a),
    b_(b),
    wf_(wf),
    wa_(wa),
    alpha_(sqrt(1/wf_ + 1/wa_)),
    beta_(sqr(cbrt(a_) + cbrt(b_)))
{}


Foam::Function2s::APIdiffCoef::APIdiffCoef
(
    const word& name,
    const dictionary& dict
)
:
    FieldFunction2<scalar, APIdiffCoef>(name)
{
    FixedList<scalar, 4> coeffs;
    Istream& is(dict.lookup(name));
    is >> coeffs;
    a_ = coeffs[0];
    b_ = coeffs[1];
    wf_ = coeffs[2];
    wa_ = coeffs[3];
    alpha_ = sqrt(1/wf_ + 1/wa_);
    beta_ = (sqr((cbrt(a_) + cbrt(b_))));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Function2s::APIdiffCoef::writeData(Ostream& os) const
{
    os.writeEntry("a", a_);
    os.writeEntry("b", b_);
    os.writeEntry("wf", wf_);
    os.writeEntry("wa", wa_);
}


// ************************************************************************* //
