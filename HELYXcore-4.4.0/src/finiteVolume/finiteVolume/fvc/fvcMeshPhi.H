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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2022-2025 Engys Ltd.

InNamespace
    Foam::fvc

Description
    Calculate the mesh motion flux and convert fluxes from absolute to relative
    and back.

SourceFiles
    fvcMeshPhi.C

\*---------------------------------------------------------------------------*/

#ifndef fvcMeshPhi_H
#define fvcMeshPhi_H

#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "dimensionedTypes/dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    template<class Type>
    tmp<surfaceScalarField> meshPhi
    (
        const VolField<Type>& vf
    );

    template<class Type>
    tmp<scalarField> meshPhi
    (
        const VolField<Type>& vf,
        const label patchi
    );

    template<class Type>
    tmp<surfaceScalarField> meshPhi
    (
        const dimensionedScalar& rho,
        const VolField<Type>& vf
    );

    template<class Type>
    tmp<surfaceScalarField> meshPhi
    (
        const volScalarField& rho,
        const VolField<Type>& vf
    );


    //- Make the given flux relative
    template<class Type>
    void makeRelative
    (
        surfaceScalarField& phi,
        const VolField<Type>& vf
    );

    //- Make the given flux relative
    template<class Type>
    void makeRelative
    (
        surfaceScalarField& phi,
        const dimensionedScalar& rho,
        const VolField<Type>& vf
    );

    //- Make the given flux relative
    template<class Type>
    void makeRelative
    (
        surfaceScalarField& phi,
        const volScalarField& rho,
        const VolField<Type>& vf
    );


    //- Make the given flux relative
    template<class Type>
    void makeRelative
    (
        surfaceScalarField& phi,
        const volScalarField& rho,
        const surfaceScalarField& rhof,
        const VolField<Type>& vf
    );


    //- Make the given flux absolute
    template<class Type>
    void makeAbsolute
    (
        surfaceScalarField& phi,
        const VolField<Type>& vf
    );

    //- Make the given flux absolute
    template<class Type>
    void makeAbsolute
    (
        surfaceScalarField& phi,
        const dimensionedScalar& rho,
        const VolField<Type>& vf
    );

    //- Make the given flux absolute
    template<class Type>
    void makeAbsolute
    (
        surfaceScalarField& phi,
        const volScalarField& rho,
        const VolField<Type>& vf
    );


    //- Return the given absolute flux in relative form
    template<class Type>
    tmp<surfaceScalarField> relative
    (
        const tmp<surfaceScalarField>& tphi,
        const VolField<Type>& vf
    );

    //- Return the given absolute flux in relative form
    template<class Type>
    tmp<surfaceScalarField> relative
    (
        const tmp<surfaceScalarField>& tphi,
        const volScalarField& rho,
        const VolField<Type>& vf
    );


    //- Return the given relative flux in absolute form
    template<class Type>
    tmp<surfaceScalarField> absolute
    (
        const tmp<surfaceScalarField>& tphi,
        const VolField<Type>& vf
    );

    //- Return the given relative flux in absolute form
    template<class Type>
    tmp<surfaceScalarField> absolute
    (
        const tmp<surfaceScalarField>& tphi,
        const volScalarField& rho,
        const VolField<Type>& vf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/fvc/fvcMeshPhiTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
