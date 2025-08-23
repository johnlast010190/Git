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
    (c) 2022 Engys Ltd

InNamespace
    Foam::fvc

Description
    Calculate the first temporal derivative.

SourceFiles
    fvcDdt.C

\*---------------------------------------------------------------------------*/


#ifndef fvcDdt_H
#define fvcDdt_H

#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "dimensionedTypes/dimensionedTypes.H"
#include "primitives/one/one.H"
#include "fields/GeometricFields/geometricZeroField/geometricZeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    template<class Type>
    tmp<VolField<Type>> ddt(const dimensioned<Type>, const fvMesh&);

    template<class Type>
    tmp<VolField<Type>> ddt(const VolField<Type>&);

    template<class Type>
    tmp<VolField<Type>> ddt(const VolField<Type>&, const word&);

    template<class Type>
    tmp<VolField<Type>> ddt(const dimensionedScalar&, const VolField<Type>&);

    template<class Type>
    tmp<VolField<Type>> ddt(const volScalarField&, const VolField<Type>&);

    template<class Type>
    tmp<VolField<Type>> ddt
    (
        const volScalarField&,
        const VolField<Type>&,
        const word&
    );

    template<class Type>
    tmp<VolField<Type>> ddt
    (
        const volScalarField&,
        const volScalarField&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> ddt
    (
        const volScalarField&,
        const volScalarField&,
        const VolField<Type>&,
        const word&
    );

    template<class Type>
    tmp<SurfaceField<Type>> ddt(const SurfaceField<Type>&);

    template<class Type>
    tmp<VolField<Type>> ddt(const one&, const VolField<Type>&);

    template<class Type>
    tmp<VolField<Type>> ddt(const VolField<Type>&, const one&);

    inline geometricZeroField ddt(const one&, const one&)
    {
        return geometricZeroField();
    }

    template<class Type>
    tmp<SurfaceField<typename Foam::flux<Type>::type>>
    ddtCorr
    (
        const VolField<Type>& U,
        const SurfaceField<Type>& Uf,
        const word interpolationName = word::null
    );

    template<class Type>
    tmp<SurfaceField<typename Foam::flux<Type>::type>>
    ddtCorr
    (
        const VolField<Type>& U,
        const SurfaceField<typename Foam::flux<Type>::type>& phi,
        const word interpolationName = word::null
    );

    template<class Type>
    tmp<SurfaceField<typename Foam::flux<Type>::type>>
    ddtCorr
    (
        const volScalarField& rho,
        const VolField<Type>& U,
        const SurfaceField<Type>& Uf,
        const word interpolationName = word::null
    );

    template<class Type>
    tmp<SurfaceField<typename Foam::flux<Type>::type>>
    ddtCorr
    (
        const volScalarField& rho,
        const VolField<Type>& U,
        const SurfaceField<typename Foam::flux<Type>::type>& phi,
        const word interpolationName = word::null
    );

    //- Convenience function for handeling two options based on Uf availability
    template<class Type>
    tmp<SurfaceField<typename Foam::flux<Type>::type>> ddtCorr
    (
        const VolField<Type>& U,
        const SurfaceField<typename Foam::flux<Type>::type>& phi,
        const autoPtr<SurfaceField<Type>>& Uf,
        const word interpolationName = word::null
    );

    //- Convenience function for handeling two options based on Uf availability
    template<class Type>
    tmp<SurfaceField<typename Foam::flux<Type>::type>> ddtCorr
    (
        const volScalarField& rho,
        const VolField<Type>& U,
        const SurfaceField<typename Foam::flux<Type>::type>& phi,
        const autoPtr<SurfaceField<Type>>& rhoUf,
        const word interpolationName = word::null
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/fvc/fvcDdt.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
