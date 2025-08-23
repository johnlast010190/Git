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
    (c) held by original author

\*---------------------------------------------------------------------------*/

#include "porosityZone.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::porosityZone::modifyDdt(fvMatrix<Type>& m) const
{
    if (porosity_ < 1)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneID_];

        scalar coeff( (1.0 + addedMassCoeff_ * ( 1 - porosity_ ) / porosity_) / porosity_ );

        forAll(cells, i)
        {
            m.diag()[cells[i]]   *= coeff;
            m.source()[cells[i]] *= coeff;
        }
    }
}

template<class RhoFieldType>
void Foam::porosityZone::addViscousInertialResistance
(
    scalarField& Udiag,
    vectorField& Usource,
    const labelList& cells,
    const scalarField& V,
    const RhoFieldType& rho,
    const scalarField& mu,
    const vectorField& U
) const
{
    const tensor& D = D_.value();
    const tensor& F = F_.value();

    forAll(cells, i)
    {
        tensor dragCoeff = mu[cells[i]]*D + (rho[cells[i]]*mag(U[cells[i]]))*F;
        scalar isoDragCoeff = tr(dragCoeff);

        Udiag[cells[i]] += V[cells[i]]*isoDragCoeff;
        Usource[cells[i]] -=
            V[cells[i]]*((dragCoeff - I*isoDragCoeff) & U[cells[i]]);
    }
}

template<class RhoFieldType>
void Foam::porosityZone::addViscousInertialResistance
(
    tensorField& AU,
    const labelList& cells,
    const RhoFieldType& rho,
    const scalarField& mu,
    const vectorField& U
) const
{
    const tensor& D = D_.value();
    const tensor& F = F_.value();

    forAll(cells, i)
    {
        AU[cells[i]] += mu[cells[i]]*D + (rho[cells[i]]*mag(U[cells[i]]))*F;
    }
}


// ************************************************************************* //
