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
    (c) 2011-2013 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "panaRotorDiskSource.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fv::panaRotorDiskSource::writeField
(
    const word& name,
    const List<Type>& values,
    const bool writeNow
) const
{
    typedef VolField<Type> fieldType;

    if (mesh_.time().outputTime() || writeNow)
    {
        tmp<fieldType> tfld
        (
            new fieldType
            (
                IOobject
                (
                    name,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensioned<Type>("zero", dimless, pTraits<Type>::zero)
            )
        );

        Field<Type>& fld = tfld.ref().ref();

        if (cells_.size() != values.size())
        {
            FatalErrorInFunction
                << "cells_.size() != values_.size()"
                << abort(FatalError);
        }

        forAll(cells_, i)
        {
            const label cellI = cells_[i];
            fld[cellI] = values[i];
        }

        tfld().write();
    }
}


// ************************************************************************* //
