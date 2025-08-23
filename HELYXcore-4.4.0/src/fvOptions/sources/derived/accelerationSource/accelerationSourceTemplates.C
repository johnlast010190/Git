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
    (c) 2018 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class AlphaRhoFieldType>
void Foam::fv::accelerationSource::add
(
    const AlphaRhoFieldType& alphaRho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    const DimensionedField<scalar, volMesh>& V = mesh_.V();

    vector a(Zero);
    const scalar t = mesh_.time().value();
    if (velOrAccel_->name() == "velocity")
    {
        const scalar dt = mesh_.time().deltaTValue();
        const vector dU = velOrAccel_->value(t) - velOrAccel_->value(t - dt);
        a = dU/mesh_.time().deltaTValue();
    }
    else if (velOrAccel_->name() == "acceleration")
    {
        a = velOrAccel_->value(t);
    }

    forAll(cells_, i)
    {
        const label c = cells_[i];
        eqn.source()[c] -= V[c]*alphaRho[c]*a;
    }
}


template<class AlphaRhoFieldType>
void Foam::fv::accelerationSource::add
(
    const AlphaRhoFieldType& alphaRho,
    fvBlockMatrix<vector>& eqn,
    const label fieldi
)
{
    const DimensionedField<scalar, volMesh>& V = mesh_.V();

    vector a(Zero);
    const scalar t = mesh_.time().value();
    if (velOrAccel_->name() == "velocity")
    {
        const scalar dt = mesh_.time().deltaTValue();
        const vector dU = velOrAccel_->value(t) - velOrAccel_->value(t - dt);
        a = dU/mesh_.time().deltaTValue();
    }
    else if (velOrAccel_->name() == "acceleration")
    {
        a = velOrAccel_->value(t);
    }

    forAll(cells_, i)
    {
        const label c = cells_[i];
        eqn.source()[c] -= V[c]*alphaRho[c]*a;
    }
}


// ************************************************************************* //
