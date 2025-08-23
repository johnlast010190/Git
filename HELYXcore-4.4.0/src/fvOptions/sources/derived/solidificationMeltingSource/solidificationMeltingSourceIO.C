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
    (c) 2014-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "solidificationMeltingSource.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fv::solidificationMeltingSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        Tmelt_ = coeffs_.lookup<scalar>("Tmelt");
        L_ = coeffs_.lookup<scalar>("L");

        coeffs_.readIfPresent("relax", relax_);

        mode_ = thermoModeTypeNames_.read(coeffs_.lookup("thermoMode"));

        rhoRef_ = coeffs_.lookup<scalar>("rhoRef");
        coeffs_.readIfPresent("T", TName_);
        coeffs_.readIfPresent("U", UName_);

        coeffs_.readIfPresent("Cu", Cu_);
        coeffs_.readIfPresent("q", q_);

        coeffs_.readIfPresent("beta", beta_);

        return true;
    }
    else
    {
        return false;
    }

    return false;
}


// ************************************************************************* //
