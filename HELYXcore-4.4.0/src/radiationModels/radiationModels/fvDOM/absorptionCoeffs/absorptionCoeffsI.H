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
    (c) 2011-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

inline bool Foam::radiationModels::absorptionCoeffs::invTemp() const
{
    return  invTemp_;
}


inline Foam::scalar Foam::radiationModels::absorptionCoeffs::Tcommon() const
{
    return  Tcommon_;
}


inline Foam::scalar Foam::radiationModels::absorptionCoeffs::Tlow() const
{
    return  Tlow_;
}


inline Foam::scalar Foam::radiationModels::absorptionCoeffs::Thigh() const
{
    return  Thigh_;
}


inline const Foam::radiationModels::absorptionCoeffs::coeffArray&
Foam::radiationModels::absorptionCoeffs::highACoeffs() const
{
    return  highACoeffs_;
}


inline const Foam::radiationModels::absorptionCoeffs::coeffArray&
Foam::radiationModels::absorptionCoeffs::lowACoeffs() const
{
    return  lowACoeffs_;
}


// ************************************************************************* //
