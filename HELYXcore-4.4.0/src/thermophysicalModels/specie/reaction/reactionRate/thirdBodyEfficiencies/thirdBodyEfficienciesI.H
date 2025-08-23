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

#include "primitives/Tuple2/Tuple2.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::thirdBodyEfficiencies::thirdBodyEfficiencies
(
    const speciesTable& species,
    const scalarList& efficiencies
)
:
    scalarField(efficiencies),
    species_(species)
{
    if (size() != species_.size())
    {
        FatalErrorInFunction
            << "number of efficiencies = " << size()
            << " is not equal to the number of species " << species_.size()
            << exit(FatalError);
    }
}


inline Foam::thirdBodyEfficiencies::thirdBodyEfficiencies
(
    const speciesTable& species,
    const dictionary& dict
)
:
    scalarField(species.size()),
    species_(species)
{
    if (dict.found("coeffs"))
    {
        List<Tuple2<word, scalar>> coeffs(dict.lookup("coeffs"));
        if (coeffs.size() != species_.size())
        {
            FatalIOErrorInFunction(dict)
                << "number of efficiencies = " << coeffs.size()
                << " is not equal to the number of species " << species_.size()
                << exit(FatalIOError);
        }

        forAll(coeffs, i)
        {
            operator[](species[coeffs[i].first()]) = coeffs[i].second();
        }
    }
    else
    {
        scalar defaultEff = dict.lookup<scalar>("defaultEfficiency");
        scalarList::operator=(defaultEff);
    }
}


// * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::thirdBodyEfficiencies::M(const scalarField& c) const
{
    scalar M = 0;
    forAll(*this, i)
    {
        M += operator[](i)*c[i];
    }

    return M;
}


inline Foam::tmp<Foam::scalarField>
Foam::thirdBodyEfficiencies::dMdc(const scalarField& c) const
{
    return *this;
}


inline void Foam::thirdBodyEfficiencies::write(Ostream& os) const
{
    List<Tuple2<word, scalar>> coeffs(species_.size());
    forAll(coeffs, i)
    {
        coeffs[i].first() = species_[i];
        coeffs[i].second() = operator[](i);
    }

    os.writeEntry("coeffs", coeffs);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const thirdBodyEfficiencies& tbes
)
{
    tbes.write(os);
    return os;
}


// ************************************************************************* //
