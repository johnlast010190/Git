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
    (c) 2010-2024 Engys Ltd.
    (c) 2012-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "cfdTools/general/porosityModel/fixedCoeff/fixedCoeff.H"
#include "fvMatrices/fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(fixedCoeff, 0);
        addToRunTimeSelectionTable(porosityModel, fixedCoeff, mesh);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::porosityModels::fixedCoeff::apply
(
    scalarField& Udiag,
    vectorField& Usource,
    const scalarField& V,
    const vectorField& U,
    const scalar rho
) const
{
    tmp<vectorField> tUrel = calculateUrel(U);
    const vectorField& Urel = tUrel();

    forAll(cellZoneIDs_, zoneI)
    {
        const tensorField& alphaZones = alpha_[zoneI];
        const tensorField& betaZones = beta_[zoneI];

        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label celli = cells[i];
            const label j = fieldIndex(i);
            const tensor Cd
                = rho*(alphaZones[j] + betaZones[j]*mag(Urel[celli]));
            const scalar isoCd = tr(Cd);

            Udiag[celli] += V[celli]*isoCd;
            Usource[celli] -= V[celli]*((Cd & Urel[celli])
                                      - ((I*isoCd) & U[celli]));
        }
    }
}


void Foam::porosityModels::fixedCoeff::apply
(
    tensorField& AU,
    vectorField& source,
    const vectorField& U,
    const scalar rho
) const
{
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

    const scalarField& V = mesh_.V();

    forAll(cellZoneIDs_, zoneI)
    {
        const tensorField& alphaZones = alpha_[zoneI];
        const tensorField& betaZones = beta_[zoneI];

        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label celli = cells[i];
            const label j = fieldIndex(i);
            const tensor alpha = alphaZones[j];
            const tensor beta = betaZones[j];

            tensor  Cd = rho*((alpha + beta*mag(Urel[celli]))).T();
            AU[celli] += Cd*V[celli];
            if (tframeU.valid())
            {
                source[celli] += Cd&tframeU()[celli]*V[celli];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::fixedCoeff::fixedCoeff
(
    const word& name,
    const word& modelType,
    const objectRegistry& obr,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    porosityModel(name, modelType, obr, mesh, dict, cellZoneName),
    alphaXYZ_("alpha", dimless/dimTime, coeffs_),
    betaXYZ_("beta", dimless/dimLength, coeffs_),
    alpha_(cellZoneIDs_.size()),
    beta_(cellZoneIDs_.size())
{
    adjustNegativeResistance(alphaXYZ_);
    adjustNegativeResistance(betaXYZ_);

    calcTransformModelData();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModels::fixedCoeff::~fixedCoeff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModels::fixedCoeff::calcTransformModelData()
{
    // The alpha coefficient as a tensor
    tensor alphaCoeff(Zero);
    alphaCoeff.xx() = alphaXYZ_.value().x();
    alphaCoeff.yy() = alphaXYZ_.value().y();
    alphaCoeff.zz() = alphaXYZ_.value().z();

    // The beta coefficient as a tensor
    tensor betaCoeff(Zero);
    betaCoeff.xx() = betaXYZ_.value().x();
    betaCoeff.yy() = betaXYZ_.value().y();
    betaCoeff.zz() = betaXYZ_.value().z();

    if (csys().uniform())
    {
        forAll(cellZoneIDs_, zonei)
        {
            alpha_[zonei].resize(1);
            beta_[zonei].resize(1);

            alpha_[zonei] = csys().transform(alphaCoeff);
            beta_[zonei] = csys().transform(betaCoeff);
        }
    }
    else
    {
        forAll(cellZoneIDs_, zonei)
        {
            const pointUIndList cc
            (
                mesh_.cellCentres(),
                mesh_.cellZones()[cellZoneIDs_[zonei]]
            );

            alpha_[zonei] = csys().transform(cc, alphaCoeff);
            beta_[zonei] = csys().transform(cc, betaCoeff);
        }
    }
}


void Foam::porosityModels::fixedCoeff::calcForce
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    vectorField& force
) const
{
    scalarField Udiag(U.size(), 0.0);
    vectorField Usource(U.size(), Zero);
    const scalarField& V = mesh_.V();
    scalar rhoRef = coeffs_.lookup<scalar>("rhoRef");

    apply(Udiag, Usource, V, U, rhoRef);

    force = Udiag*U - Usource;
}


void Foam::porosityModels::fixedCoeff::correct
(
    fvVectorMatrix& UEqn
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    scalar rho = 1.0;
    if (UEqn.dimensions() == dimForce)
    {
        rho = coeffs_.lookup<scalar>("rhoRef");
    }

    apply(Udiag, Usource, V, U, rho);
}


void Foam::porosityModels::fixedCoeff::correct
(
    fvBlockMatrix<vector>& UEqn
) const
{
    const volVectorField& U = UEqn.psi();
    vectorField& Usource = UEqn.source();
    typename CoeffField<vector>::squareTypeField& blockDiag =
        UEqn.diag().asSquare();

    scalar rho = 1.0;
    if (UEqn.dimensionSets()[0] == dimForce)
    {
        rho = coeffs_.lookup<scalar>("rhoRef");
    }

    apply(blockDiag, Usource, U, rho);
}


void Foam::porosityModels::fixedCoeff::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField&,
    const volScalarField&
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    scalar rho = 1.0;
    if (UEqn.dimensions() == dimForce)
    {
        rho = coeffs_.lookup<scalar>("rhoRef");
    }

    apply(Udiag, Usource, V, U, rho);
}


void Foam::porosityModels::fixedCoeff::correct
(
    fvBlockMatrix<vector>& UEqn,
    const volScalarField&,
    const volScalarField&
) const
{
    const volVectorField& U = UEqn.psi();

    typename CoeffField<vector>::squareTypeField& blockDiag =
        UEqn.diag().asSquare();
    vectorField& Usource = UEqn.source();
    scalar rho = 1.0;
    if (UEqn.dimensionSets()[0] == dimForce)
    {
        rho = coeffs_.lookup<scalar>("rhoRef");
    }

    apply(blockDiag, Usource, U, rho);
}


void Foam::porosityModels::fixedCoeff::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const vectorField& U = UEqn.psi();
    vectorField nullSource (AU.size(), 0);
    scalar rho = 1.0;
    if (UEqn.dimensions() == dimForce)
    {
        rho = coeffs_.lookup<scalar>("rhoRef");
    }

    apply(AU, nullSource, U, rho);
}


void Foam::porosityModels::fixedCoeff::adjointCorrect
(
    fvVectorMatrix& UaEqn,
    const volVectorField& Uprimal
) const
{}


void Foam::porosityModels::fixedCoeff::adjointCorrect
(
    fvVectorMatrix& UaEqn,
    const volScalarField& rho,
    const volScalarField& mu,
    const volVectorField& Uprimal
) const
{}


void Foam::porosityModels::fixedCoeff::adjointCorrect
(
    fvBlockMatrix<vector>&,
    const volVectorField&
) const
{
   NotImplemented;
}


void Foam::porosityModels::fixedCoeff::adjointCorrect
(
    fvBlockMatrix<vector>&,
    const volScalarField&,
    const volScalarField&,
    const volVectorField&
) const
{
   NotImplemented;
}


bool Foam::porosityModels::fixedCoeff::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);

    return true;
}


// ************************************************************************* //
