/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) ICE Stroemungsfoschungs GmbH

Class
    Foam::fv::matrixChangeBaseFvOption

Contributors/Copyright:
    2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef matrixChangeBaseFvOptionI_H
#define matrixChangeBaseFvOptionI_H

#include "include/swak.H"

namespace Foam
{
namespace fv
{

template<class T>
tmp<VolField<T>>
matrixChangeBaseFvOption::calcResiduum(const fvMatrix<T> &matrix) const
{
    tmp<VolField<T>> tres(
        new VolField<T>(
            IOobject(
                "tmp"+this->fieldName(),
                matrix.psi().mesh().time().timeName(),
                matrix.psi().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            matrix.psi().mesh(),
            dimensioned<T>(
                "no",
                matrix.dimensions()/dimVolume,
                pTraits<T>::zero
            )
        )
    );
    VolField<T> &res=
        const_cast<VolField<T>&>(tres());

#ifdef FOAM_NO_DIMENSIONEDINTERNAL_IN_GEOMETRIC
    const_cast<DimensionedField<T,volMesh>&>(res.internalField())=
#else
    res.internalField() =
#endif
        (matrix & matrix.psi());
#ifdef FOAM_NO_DIMENSIONEDINTERNAL_IN_GEOMETRIC
    const_cast<Field<T>&>(res.internalField().field())-=
#else
    res.internalField() -=
#endif
        matrix.source();

    return tres;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
