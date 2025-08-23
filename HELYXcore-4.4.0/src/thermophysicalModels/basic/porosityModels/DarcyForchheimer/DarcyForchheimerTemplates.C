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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2010-2018 Engys Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::porosityModels::DarcyForchheimer::apply
(
    scalarField& Udiag,
    vectorField& Usource,
    const scalarField& V,
    const RhoFieldType& rho,
    const scalarField& mu,
    const vectorField& U
) const
{
    tmp<vectorField> tUrel = calculateUrel(U);
    const vectorField& Urel = tUrel();

    forAll(cellZoneIDs_, zoneI)
    {
        const tensorField& dZones = D_[zoneI];
        const tensorField& fZones = F_[zoneI];

        const scalarField& ddnZones = Dcorr_[zoneI];
        const scalarField& ffnZones = Fcorr_[zoneI];

        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label celli = cells[i];
            const label j = fieldIndex(i);

            tensor Cd =
                    mu[celli]*dZones[j]*ddnZones[i] +
                    (
                        rho[celli] *
                        mag(Urel[celli])
                    ) *
                    fZones[j] *
                    ffnZones[i];

            if (!defaultMode_)
            {
                Cd.xx() = max(Cd.xx(), 0.0);
                Cd.yy() = max(Cd.yy(), 0.0);
                Cd.zz() = max(Cd.zz(), 0.0);
            }

            const scalar isoCd = tr(Cd);

            Udiag[celli] += V[celli]*isoCd;
            Usource[celli] -= V[celli]*((Cd & Urel[celli])
                                      - ((I*isoCd) & U[celli]));
        }
    }
}


template<class RhoFieldType>
void Foam::porosityModels::DarcyForchheimer::apply
(
    tensorField& AU,
    vectorField& source,
    const RhoFieldType& rho,
    const scalarField& mu,
    const vectorField& U
) const
{
    const scalarField& V = mesh_.V();

    tmp<vectorField> tUrel(new vectorField(U));
    vectorField& Urel = tUrel.ref();
    substractMeshVelocity(Urel);

    tmp<vectorField> tframeU;
    if (coorFramePtr_)
    {
        tframeU = new vectorField(U.size(), Zero);
        vectorField& frameU = tframeU.ref();
        frameVelocity(frameU, true);
        Urel -= tframeU();
    }


    forAll(cellZoneIDs_, zoneI)
    {
        const tensorField& dZones = D_[zoneI];
        const tensorField& fZones = F_[zoneI];

        const scalarField& ddnZones = Dcorr_[zoneI];
        const scalarField& ffnZones = Fcorr_[zoneI];

        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label cI = cells[i];
            const label j = fieldIndex(i);
            const tensor D = dZones[j]*ddnZones[i];
            const tensor F = fZones[j]*ffnZones[i];

            tensor Cd = mu[cI]*D + (rho[cI]*mag(Urel[cI]))*F;

            if (!defaultMode_)
            {
                Cd.xx() = max(Cd.xx(), 0.0);
                Cd.yy() = max(Cd.yy(), 0.0);
                Cd.zz() = max(Cd.zz(), 0.0);
            }

            AU[cI] += Cd*V[cI];

            if (coorFramePtr_)
            {
                source[cI] += Cd&tframeU()[cI]*V[cI];
            }
        }
    }
}

template<class RhoFieldType>
void Foam::porosityModels::DarcyForchheimer::adjointApply
(
    scalarField& Udiag,
    vectorField& Usource,
    const scalarField& V,
    const RhoFieldType& rho,
    const scalarField& mu,
    const vectorField& Uprimal,
    const vectorField& U
) const
{
    vectorField UpRel = Uprimal;

    frameVelocity(UpRel, false);

    forAll(cellZoneIDs_, zoneI)
    {
        tensorField dZonest = D_[zoneI].T()();
        tensorField fZonest = F_[zoneI].T()();
        const tensorField& fZones = F_[zoneI];

        const scalarField& ddnZones = Dcorr_[zoneI];
        const scalarField& ffnZones = Fcorr_[zoneI];

        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label cellI = cells[i];
            const label j = fieldIndex(i);

            const vector& Upr(UpRel[cellI]);
            scalar magUprimal = max(mag(Upr), SMALL);

            tensor dragCoeff =
                    mu[cellI]*dZonest[j]*ddnZones[i]
                  + (rho[cellI]*magUprimal)*fZonest[j]*ffnZones[i]
                  + (rho[cellI]/(magUprimal)
                  * (fZones[j]*ffnZones[i] & sqr(Upr)).T());

            if (!defaultMode_)
            {
                dragCoeff.xx() = max(dragCoeff.xx(), 0.0);
                dragCoeff.yy() = max(dragCoeff.yy(), 0.0);
                dragCoeff.zz() = max(dragCoeff.zz(), 0.0);
            }

            const scalar isoDragCoeff = tr(dragCoeff);

            Udiag[cellI] += V[cellI]*isoDragCoeff;
            Usource[cellI] -=
                (V[cellI]*((dragCoeff - I*isoDragCoeff) & U[cellI]));
        }
    }
}

// ************************************************************************* //
