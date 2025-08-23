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
    (c) 2010-2022 Engys Ltd

\*---------------------------------------------------------------------------*/

inline const Foam::volScalarField&
Foam::radiationModels::radiativeIntensityRay::I() const
{
    return I_;
}

inline const Foam::volScalarField::Boundary&
Foam::radiationModels::radiativeIntensityRay::qr() const
{
    return qr_();
}


inline Foam::volScalarField::Boundary& Foam::radiationModels::radiativeIntensityRay::qr()
{
    return qr_();
}

inline const Foam::volScalarField::Boundary& Foam::radiationModels::
radiativeIntensityRay::qin() const
{
    return qin_();
}


inline Foam::volScalarField::Boundary& Foam::radiationModels::radiativeIntensityRay::qin()
{
    return qin_();
}


inline const Foam::volScalarField::Boundary& Foam::radiationModels::
radiativeIntensityRay::qem() const
{
    return qem_();
}


inline Foam::volScalarField::Boundary& Foam::radiationModels::radiativeIntensityRay::qem()
{
    return qem_();
}

inline const Foam::volScalarField::Boundary& Foam::radiationModels::
radiativeIntensityRay::qg() const
{
    return qg_();
}


inline Foam::volScalarField::Boundary& Foam::radiationModels::radiativeIntensityRay::qg()
{
    return qg_();
}


inline const Foam::vector& Foam::radiationModels::radiativeIntensityRay::d() const
{
    return d_;
}


inline Foam::vector& Foam::radiationModels::radiativeIntensityRay::d()
{
    return d_;
}


inline const Foam::vector& Foam::radiationModels::radiativeIntensityRay::dAve() const
{
    return dAve_;
}


inline Foam::vector& Foam::radiationModels::radiativeIntensityRay::dAve()
{
    return dAve_;
}


inline Foam::scalar Foam::radiationModels::radiativeIntensityRay::nLambda() const
{
    return nLambda_;
}


inline Foam::scalar Foam::radiationModels::radiativeIntensityRay::phi() const
{
    return phi_;
}


inline Foam::scalar Foam::radiationModels::radiativeIntensityRay::theta() const
{
    return theta_;
}


inline Foam::scalar Foam::radiationModels::radiativeIntensityRay::omega() const
{
    return omega_;
}

inline Foam::word Foam::radiationModels::radiativeIntensityRay::rayName() const
{
    return word("I" + name(myRayId_));
}

inline const Foam::volScalarField&
Foam::radiationModels::radiativeIntensityRay::ILambda
(
    const label lambdaI
) const
{
    return ILambda_[lambdaI];
}


// ************************************************************************* //
