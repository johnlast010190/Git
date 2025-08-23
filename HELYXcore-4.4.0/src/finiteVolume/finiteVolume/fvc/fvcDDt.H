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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2020 Engys Ltd.

InNamespace
    Foam::fvc

Description
    Calculate the substantive (total) derivative.

Usage
    - DDt(phi,psi): calculates D(psi)/Dt = d(psi)/dt + u*div(psi)
    - DDt(rho,phi,psi): calculates D(psi)/Dt = d(psi)/dt) + div(phi*u)

SourceFiles
    fvcDDt.C

\*---------------------------------------------------------------------------*/


#ifndef fvcDDt_H
#define fvcDDt_H

#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    template<class Type>
    tmp<VolField<Type>> DDt
    (
        const surfaceScalarField&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> DDt
    (
        const geometricOneField& rho,
        const surfaceScalarField& phi,
        const VolField<Type>& psi
    );

    template<class Type>
    tmp<VolField<Type>> DDt
    (
        const volScalarField& rho,
        const surfaceScalarField& phi,
        const VolField<Type>& psi
    );

    template<class Type>
    tmp<VolField<Type>> DDt
    (
        const tmp<surfaceScalarField>&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> DDt
    (
        const geometricOneField& rho,
        const tmp<surfaceScalarField>&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> DDt
    (
        const volScalarField& rho,
        const tmp<surfaceScalarField>&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> DDt
    (
        const tmp<geometricOneField>&,
        const surfaceScalarField& phi,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> DDt
    (
        const tmp<volScalarField>&,
        const surfaceScalarField& phi,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> DDt
    (
        const tmp<geometricOneField>&,
        const tmp<surfaceScalarField>&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> DDt
    (
        const tmp<volScalarField>&,
        const tmp<surfaceScalarField>&,
        const VolField<Type>&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/fvc/fvcDDt.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
