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
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SolidThermo, class GasThermo>
inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::pyrolysisChemistryModel<SolidThermo, GasThermo>::RRg()
{
    return RRg_;
}


template<class SolidThermo, class GasThermo>
inline Foam::label
Foam::pyrolysisChemistryModel<SolidThermo, GasThermo>::nSpecie() const
{
    return nSpecie_;
}


template<class SolidThermo, class GasThermo>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::pyrolysisChemistryModel<SolidThermo, GasThermo>::RRg(const label i) const
{
    return RRg_[i];
}


template<class SolidThermo, class GasThermo>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::pyrolysisChemistryModel<SolidThermo, GasThermo>::RRg() const
{
    tmp<volScalarField::Internal> tRRg
    (
        volScalarField::Internal::New
        (
            "RRg",
            this->mesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, 0)
        )
    );

    if (this->chemistry_)
    {
        volScalarField::Internal& RRg = tRRg.ref();
        for (label i=0; i < nGases_; i++)
        {
            RRg += RRg_[i];
        }
    }

    return tRRg;
}


// ************************************************************************* //
