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
    (c) 2016 OpenCFD Ltd.
    (c) 2013-2019 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SolidThermo>
inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::solidChemistryModel<SolidThermo>::RRs()
{
    return RRs_;
}


template<class SolidThermo>
inline const Foam::ReactionList<SolidThermo>&
Foam::solidChemistryModel<SolidThermo>::reactions() const
{
    return reactions_;
}


template<class SolidThermo>
inline Foam::label
Foam::solidChemistryModel<SolidThermo>::nReaction() const
{
    return nReaction_;
}


template<class SolidThermo>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::solidChemistryModel<SolidThermo>::RRs
(
    const label i
) const
{
    return RRs_[i];
}


template<class SolidThermo>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::solidChemistryModel<SolidThermo>::RRs() const
{
    tmp<volScalarField::Internal> tRRs
    (
        volScalarField::Internal::New
        (
            "RRs",
            this->mesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, 0)
        )
    );

    if (this->chemistry_)
    {
        volScalarField::Internal& RRs = tRRs.ref();
        for (label i=0; i < nSolids_; i++)
        {
            RRs += RRs_[i];
        }
    }
    return tRRs;
}


template<class SolidThermo>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::solidChemistryModel<SolidThermo>::RRsHs() const
{
    tmp<DimensionedField<scalar, volMesh>> tRRsHs
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "RRsHs",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
        )
    );

    if (this->chemistry_)
    {
        DimensionedField<scalar, volMesh>& RRs = tRRsHs.ref();
        if (basicThermo::dictName == basicThermo::matDictName)
        {
            for (label i=0; i < nSolids_; i++)
            {
                RRs.field() += RRs_[i].field()*hsModels_[i].primitiveField();
            }
        }
        else
        {
            const volScalarField& T = this->solidThermo().T();
            const volScalarField& p = this->solidThermo().p();
            for (label i=0; i < nSolids_; i++)
            {
                forAll(RRs, celli)
                {
                    RRs[celli] +=
                        RRs_[i][celli]*solidThermo_[i].hs(p[celli], T[celli]);
                }
            }
        }
    }
    return tRRsHs;
}


// ************************************************************************* //
