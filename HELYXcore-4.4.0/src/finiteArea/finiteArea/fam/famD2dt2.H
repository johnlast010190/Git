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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Engys Ltd.

Namespace
    fam

Description
    Calculate the matrix for the second temporal derivative.

SourceFiles
    famD2dt2.C

\*---------------------------------------------------------------------------*/

#ifndef famD2dt2_H
#define famD2dt2_H

#include "fields/areaFields/areaFieldsFwd.H"
#include "faMatrices/faMatrix/faMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fam functions Declaration
\*---------------------------------------------------------------------------*/

namespace fam
{
    template<class Type>
    tmp<faMatrix<Type>> d2dt2
    (
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<faMatrix<Type>> d2dt2
    (
        const dimensionedScalar&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<faMatrix<Type>> d2dt2
    (
        const areaScalarField&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "finiteArea/fam/famD2dt2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
