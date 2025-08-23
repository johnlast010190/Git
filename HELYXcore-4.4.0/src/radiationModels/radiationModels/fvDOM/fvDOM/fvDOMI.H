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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2010-2020 Engys Ltd

\*---------------------------------------------------------------------------*/

inline const Foam::radiationModels::radiativeIntensityRay&
Foam::radiationModels::fvDOM::IRay(const label rayI) const
{
    return  IRay_[rayI];
}


inline const Foam::volScalarField&
Foam::radiationModels::fvDOM::IRayLambda
(
    const label rayI,
    const label lambdaI
) const
{
    return IRay_[rayI].ILambda(lambdaI);
}


inline Foam::label Foam::radiationModels::fvDOM::nTheta() const
{
    return nTheta_;
}


inline Foam::label Foam::radiationModels::fvDOM::nPhi() const
{
    return nPhi_;
}


inline Foam::label Foam::radiationModels::fvDOM::nRay() const
{
    return nRay_;
}


inline Foam::label Foam::radiationModels::fvDOM::nLambda() const
{
    return nLambda_;
}


inline Foam::label Foam::radiationModels::fvDOM::nBands() const
{
    return nLambda_;
}


inline const Foam::volScalarField& Foam::radiationModels::fvDOM::a() const
{
    return a_;
}


inline const Foam::volScalarField& Foam::radiationModels::fvDOM::aLambda
(
    const label lambdaI
) const
{
    return aLambda_[lambdaI];
}


inline const Foam::volScalarField& Foam::radiationModels::fvDOM::G() const
{
    return G_;
}


inline const Foam::volScalarField& Foam::radiationModels::fvDOM::qr() const
{
    return qr_;
}

inline const Foam::volScalarField& Foam::radiationModels::fvDOM::qin() const
{
    return qin_;
}

//inline const Foam::volScalarField& Foam::radiationModels::fvDOM::qem() const
//{
//    return qem_;
//}

inline const Foam::volScalarField& Foam::radiationModels::fvDOM::qg() const
{
    return qg_;
}


inline const Foam::radiationModels::blackBodyEmission&
Foam::radiationModels::fvDOM::blackBody() const
{
    return blackBody_;
}


inline const Foam::waveSolverCoeffs& Foam::radiationModels::fvDOM::rayWaveCoeffs
(
    const label rayId,
    const label lambdaI
) const
{
    return rayWaveCoeffs_[lambdaI][rayId];
}

inline Foam::scalar Foam::radiationModels::fvDOM::omegaMax() const
{
    return omegaMax_;
}


inline Foam::vector Foam::radiationModels::fvDOM::meshOrientation() const
{
    return meshOrientation_;
}


inline bool Foam::radiationModels::fvDOM::useExternalBeam() const
{
    return useExternalBeam_;
}


inline const Foam::radiationModels::boundaryRadiationProperties&
Foam::radiationModels::fvDOM::boundaryProperties() const
{
    return boundaryRadiationProperties::New(mesh_, TRef_.value());
}

inline const Foam::scalarList& Foam::radiationModels::fvDOM::
spectralDistribution() const
{
    return spectralDistribution_;
}

// ************************************************************************* //
