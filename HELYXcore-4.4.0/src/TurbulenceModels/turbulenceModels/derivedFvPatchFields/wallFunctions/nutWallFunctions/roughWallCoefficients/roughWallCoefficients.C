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
    (c) 2020 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "roughWallCoefficients.H"

void Foam::roughWallCoefficients::checkOldFormat(const dictionary& dict)
{
    if (dict.found("Ks") && !roughnessIsActive())
    {
        roughnessHeight_ = dict.lookupOrDefault<scalar>("Ks", 0.);
        roughnessConstant_ = dict.lookupOrDefault<scalar>("Cs", 0.);
        WarningInFunction
            << "Old style definition of roughness coefficients is used"<<nl
            << "Ks is now substituted by roughnessHeight"<<nl
            << "and Cs by roughnessConstant."<<nl
            << "Please consider using the new format "<<endl;
    }
    else if (dict.found("z0") && !roughnessIsActive())
    {
        roughnessHeight_ = dict.lookupOrDefault<scalar>("z0", 0.);
        WarningInFunction
        << "Old style definition of roughness coefficients is used"<<nl
        << "z0 is now substituted by roughnessHeight."<<nl
        << "Please consider using the new format "<<endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::roughWallCoefficients::roughWallCoefficients(const Foam::dictionary &dict)
:
roughnessHeight_(dict.lookupOrDefault<scalar>("roughnessHeight", 0.)),
roughnessConstant_(dict.lookupOrDefault<scalar>("roughnessConstant", 0.5)),
roughnessFactor_(dict.lookupOrDefault<scalar>("roughnessFactor", 1.))
{
    checkOldFormat(dict);
}

Foam::roughWallCoefficients::roughWallCoefficients()
:
roughnessHeight_(Zero),
roughnessConstant_(Zero),
roughnessFactor_(Zero)
{}

// * * * * * * * * * * * * * * * * Member Functions   * * * * * * * * * * * *//
void Foam::roughWallCoefficients::writeLocalEntries(Foam::Ostream &os) const
{
    if (roughnessIsActive())
    {
        os.writeEntry("roughnessHeight", roughnessHeight_);
        os.writeEntry("roughnessConstant", roughnessConstant_);
        os.writeEntry("roughnessFactor", roughnessFactor_);
    }
}
bool Foam::roughWallCoefficients::roughnessIsActive() const
{
    if (roughnessHeight_ > 0.)
    {
        return true;
    }

    return false;
}