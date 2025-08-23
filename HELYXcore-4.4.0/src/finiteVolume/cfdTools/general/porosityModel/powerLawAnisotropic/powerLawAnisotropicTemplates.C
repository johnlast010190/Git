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
    (c) 2018 Engys Ltd.
    (c) 2012 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::porosityModels::powerLawAnisotropic::apply
(
    scalarField& Udiag,
    vectorField& Usource,
    const scalarField& V,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    tmp<vectorField> tUrel = calculateUrel(U);
    const vectorField& Urel = tUrel();

    const scalar  Bm1over2 = (B_ - 1.0)/2.0;

    forAll(cellZoneIDs_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label celli = cells[i];

            const tensor Cd = rho[celli]*C_*pow(magSqr(Urel[celli]), Bm1over2);

            const scalar isoCd = tr(Cd);

            Udiag[celli] += V[celli]*isoCd;
            Usource[celli] -= V[celli]*((Cd & Urel[celli])  - ((I*isoCd) & U[celli]));
        }
    }
}


template<class RhoFieldType>
void Foam::porosityModels::powerLawAnisotropic::apply
(
    tensorField& AU,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    if (validParentFrame())
    {
        FatalErrorInFunction
            << "Tensorial anisotropic powerLaw not supported in moving "
            << "reference frame." << exit(FatalError);
    }

    const scalar  Bm1over2 = (B_ - 1.0)/2.0;

    forAll(cellZoneIDs_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label celli = cells[i];
            AU[celli] += rho[celli]*C_*pow(magSqr(U[celli]), Bm1over2);
        }
    }
}


template<class RhoFieldType>
void Foam::porosityModels::powerLawAnisotropic::adjointApply
(
    scalarField& Udiag,
    vectorField& Usource,
    const scalarField& V,
    const RhoFieldType& rho,
    const vectorField& Uprimal,
    const vectorField& U
) const
{

    vectorField UpRel = Uprimal;
    frameVelocity(UpRel, false);
    const scalar  Bm1over2 = (B_ - 1.0)/2.0;
    const tensor& Ct = C_.T();

    forAll(cellZoneIDs_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];
        forAll(cells, i)
        {
            const vector& Upr(UpRel[cells[i]]);

            const tensor dragCoeff
                = rho[cells[i]]*Ct*pow(magSqr(Upr), Bm1over2);

            const scalar isoDragCoeff = tr(dragCoeff);

            Udiag[cells[i]] += V[cells[i]]*isoDragCoeff;
            Usource[cells[i]] -=
                (V[cells[i]]*((dragCoeff - I*isoDragCoeff) & U[cells[i]]));

        }
    }
}

// ************************************************************************* //
