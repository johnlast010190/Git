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
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "global/constants/thermodynamic/thermodynamicConstants.H"
#include "global/constants/physicoChemical/physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constant
{
namespace thermodynamic
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Note: the 1e3 converts from /mol to /kmol for consistency with the
    // SI choice of kg rather than g for mass.
    // This is not appropriate for USCS and will be changed to an entry in
    // the DimensionedConstants dictionary in etc/controlDict
    const scalar RR = 1e3*physicoChemical::R.value();

    const scalar Pstd = standard::Pstd.value();
    const scalar Tstd = standard::Tstd.value();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermodynamic
} // End namespace constant
} // End namespace Foam

// ************************************************************************* //
