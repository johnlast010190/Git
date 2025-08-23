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

InNamespace
    Foam::fvc

Description
    Volume integrate volField creating a volField.

    Volume integrate volField over the whole domain
    creating a dimensioned\<Type\>

SourceFiles
    fvcVolumeIntegrate.C

\*---------------------------------------------------------------------------*/


#ifndef fvcVolumeIntegrate_H
#define fvcVolumeIntegrate_H

#include "fields/volFields/volFieldsFwd.H"
#include "fields/Fields/primitiveFieldsFwd.H"
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
    tmp<Field<Type>> volumeIntegrate
    (
        const VolField<Type>&
    );

    template<class Type>
    tmp<Field<Type>> volumeIntegrate
    (
        const tmp<VolField<Type>>&
    );


    template<class Type>
    tmp<Field<Type>> volumeIntegrate
    (
        const DimensionedField<Type, volMesh>&
    );

    template<class Type>
    tmp<Field<Type>> volumeIntegrate
    (
        const tmp<DimensionedField<Type, volMesh>>&
    );


    template<class Type>
    dimensioned<Type> domainIntegrate
    (
        const VolField<Type>&
    );

    template<class Type>
    dimensioned<Type> domainIntegrate
    (
        const tmp<VolField<Type>>&
    );


    template<class Type>
    dimensioned<Type> domainIntegrate
    (
        const DimensionedField<Type, volMesh>&
    );

    template<class Type>
    dimensioned<Type> domainIntegrate
    (
        const tmp<DimensionedField<Type, volMesh>>&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/fvc/fvcVolumeIntegrate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
