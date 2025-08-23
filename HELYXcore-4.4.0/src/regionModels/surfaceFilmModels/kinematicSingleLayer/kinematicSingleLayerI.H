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
    (c) 2011-2021 OpenFOAM Foundation
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "finiteVolume/fvc/fvcSurfaceIntegrate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


inline const dimensionedScalar& kinematicSingleLayer::deltaSmall() const
{
    return deltaSmall_;
}


inline const volScalarField& kinematicSingleLayer::rho() const
{
    return thermo_->rho();
}


inline const volScalarField& kinematicSingleLayer::mu() const
{
    return mu_;
}


inline const volScalarField& kinematicSingleLayer::delta() const
{
    return delta_;
}


inline const volScalarField& kinematicSingleLayer::alpha() const
{
    return alpha_;
}


inline const volVectorField& kinematicSingleLayer::U() const
{
    return U_;
}


inline tmp<volVectorField::Internal> kinematicSingleLayer::Us() const
{
    return turbulence_->Us();
}


inline const surfaceScalarField& kinematicSingleLayer::phi() const
{
    return phi_;
}


inline const surfaceScalarField& kinematicSingleLayer::phiU() const
{
    return phiU_;
}


inline const volScalarField::Internal&
kinematicSingleLayer::continuityErr() const
{
    return continuityErr_;
}


inline const volScalarField& kinematicSingleLayer::coverage() const
{
    return coverage_;
}


inline volVectorField& kinematicSingleLayer::USpPrimary()
{
    return USpPrimary_;
}


inline volScalarField& kinematicSingleLayer::pSpPrimary()
{
    return pSpPrimary_;
}


inline volScalarField& kinematicSingleLayer::rhoSpPrimary()
{
    return rhoSpPrimary_;
}


inline volVectorField::Internal& kinematicSingleLayer::USp()
{
    return USp_;
}


inline volScalarField::Internal& kinematicSingleLayer::pSp()
{
    return pSp_;
}


inline volScalarField::Internal& kinematicSingleLayer::rhoSp()
{
    return rhoSp_;
}


inline const volVectorField::Internal& kinematicSingleLayer::USp() const
{
    return USp_;
}


inline const volScalarField::Internal& kinematicSingleLayer::pSp() const
{
    return pSp_;
}


inline const volScalarField::Internal& kinematicSingleLayer::rhoSp() const
{
    return rhoSp_;
}


inline const volVectorField& kinematicSingleLayer::UPrimary() const
{
    return UPrimary_;
}


inline const volScalarField& kinematicSingleLayer::pPrimary() const
{
    return thermo_->p();
}


inline const volScalarField& kinematicSingleLayer::rhoPrimary() const
{
    return rhoPrimary_;
}


inline const volScalarField& kinematicSingleLayer::muPrimary() const
{
    return muPrimary_;
}


inline const rhoThermo& kinematicSingleLayer::thermo() const
{
    return thermo_();
}


inline ejectionModelList& kinematicSingleLayer::ejection()
{
    return ejection_;
}


inline transferModelList& kinematicSingleLayer::transfer()
{
    return transfer_;
}


inline const turbulenceModel& kinematicSingleLayer::turbulence() const
{
    return turbulence_();
}


inline tmp<volScalarField::Internal> kinematicSingleLayer::mass() const
{
    return rho()()*delta_()*magSf();
}


inline tmp<volScalarField::Internal> kinematicSingleLayer::deltaMass() const
{
    return rhoSp_*magSf()*time().deltaT();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
