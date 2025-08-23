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
    (c) 2021 Engys Ltd.
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "materialsSpecie.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::materialsSpecie>
Foam::materialsSpecie::clone() const
{
    return autoPtr<materialsSpecie>(new materialsSpecie(*this));
}


inline Foam::autoPtr<Foam::materialsSpecie> Foam::materialsSpecie::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& materialsSpecieName,
    const word& name
)
{
    return autoPtr<materialsSpecie>
    (
        new materialsSpecie(obr, dict, phaseName, materialsSpecieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::materialsSpecie::W() const
{
    return molWeight_;
}


inline Foam::scalar Foam::materialsSpecie::Y() const
{
    return Y_;
}


inline Foam::scalar Foam::materialsSpecie::R() const
{
    return R_;
}


// ************************************************************************* //
