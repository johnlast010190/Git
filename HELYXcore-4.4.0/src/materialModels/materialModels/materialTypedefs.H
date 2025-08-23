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
    (c) 2021 Engys Ltd.

\*---------------------------------------------------------------------------*/

#ifndef materialTypedefs_H
#define materialTypedefs_H

#include "materialModels/baseModels/baseModels.H"
#include "general/referenceFields/referenceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef HashTable<baseModels<scalar>*> matScalarTable;
    typedef HashTable<baseModels<vector>*> matVectorTable;
    typedef HashTable<baseModels<tensor>*> matTensorTable;

    typedef HashTable<matScalarTable> matScalarTableTable;
    typedef HashTable<matVectorTable> matVectorTableTable;
    typedef HashTable<matTensorTable> matTensorTableTable;

    typedef UPtrList<baseModels<scalar>> matScalarUTable;
    typedef UPtrList<baseModels<vector>> matVectorUTable;
    typedef UPtrList<baseModels<tensor>> matTensorUTable;

    typedef referenceFields<scalar> refScalarField;
    typedef referenceFields<vector> refVectorField;
    typedef referenceFields<tensor> refTensorField;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
