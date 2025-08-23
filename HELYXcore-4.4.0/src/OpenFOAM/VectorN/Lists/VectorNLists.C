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
    (c) 2010 Ivor Clifford

Description
    Specialisation of List<T> for VectorN and TensorN types.

\*---------------------------------------------------------------------------*/

#include "VectorN/Lists/VectorNLists.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#define makeListType(type, Type, args...)       \
    defineCompoundTypeName(List<type>,             \
        type##List);                            \
    addCompoundToRunTimeSelectionTable(            \
        List<type>, type##List);                \
                                                \
    defineTemplateTypeNameAndDebugWithName(        \
        type##IOList, #type"List", 0);            \
                                                \
    defineTemplateTypeNameAndDebugWithName        \
    (                                            \
        type##ListIOList,                        \
        #type"ListList",                        \
        0                                        \
    );                                            \
                                                \
    defineTemplateTypeNameAndDebugWithName        \
    (                                            \
        type##ListCompactIOList,                \
        #type"ListCompactList",                    \
        0                                        \
    );                                            \

forAllVectorNTypes(makeListType)
forAllTensorNTypes(makeListType)
forAllDiagTensorNTypes(makeListType)
forAllSphericalTensorNTypes(makeListType)

#undef makeListType

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
