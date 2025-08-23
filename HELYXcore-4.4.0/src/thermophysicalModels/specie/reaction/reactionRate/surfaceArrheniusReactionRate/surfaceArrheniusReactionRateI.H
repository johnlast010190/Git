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
    (c) 2019-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::surfaceArrheniusReactionRate::surfaceArrheniusReactionRate
(
    const speciesTable& species,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    ArrheniusReactionRate(species, dict),
    AvName_(dict.lookup(wordList({"Av", "a"}))),
    obr_(obr),
    tAv_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::surfaceArrheniusReactionRate::preEvaluate() const
{
    ArrheniusReactionRate::preEvaluate();

    tAv_ = obr_.lookupObject<volScalarField::Internal>(AvName_);
}


inline void Foam::surfaceArrheniusReactionRate::postEvaluate() const
{
    ArrheniusReactionRate::postEvaluate();

    tAv_ = tmp<volScalarField::Internal>(nullptr);
}


inline Foam::scalar Foam::surfaceArrheniusReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return ArrheniusReactionRate::operator()(p, T, c, li)*tAv_()[li];
}


inline Foam::scalar Foam::surfaceArrheniusReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return ArrheniusReactionRate::ddT(p, T, c, li)*tAv_()[li];
}


inline void Foam::surfaceArrheniusReactionRate::write(Ostream& os) const
{
    ArrheniusReactionRate::write(os);
    os.writeEntry("a", AvName_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const surfaceArrheniusReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //