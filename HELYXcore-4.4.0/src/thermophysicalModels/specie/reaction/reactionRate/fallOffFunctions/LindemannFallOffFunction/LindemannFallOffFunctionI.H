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
    (c) 2011-2021 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::LindemannFallOffFunction::LindemannFallOffFunction()
{}


inline Foam::LindemannFallOffFunction::LindemannFallOffFunction
(
    const dictionary& dict
)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::LindemannFallOffFunction::operator()
(
    const scalar,
    const scalar
) const
{
    return 1.0;
}


inline Foam::scalar Foam::LindemannFallOffFunction::ddT
(
    const scalar,
    const scalar,
    const scalar
) const
{
    return 0;
}


inline Foam::scalar Foam::LindemannFallOffFunction::ddPr
(
    const scalar,
    const scalar,
    const scalar
) const
{
    return 0;
}


inline void Foam::LindemannFallOffFunction::write(Ostream& os) const
{}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::LindemannFallOffFunction&
)
{
    return os;
}


// ************************************************************************* //
