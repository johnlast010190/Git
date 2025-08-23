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
    (c) 2019-2022 OpenFOAM Foundation
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/Fields/fieldMappers/setSizeFieldMapper/setSizeFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setSizeFieldMapper::setSizeFieldMapper(const label size)
:
    size_(size)
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

FOR_ALL_FIELD_TYPES(DEFINE_FIELD_MAPPER_OPERATOR, , setSizeFieldMapper)


DEFINE_FIELD_MAPPER_OPERATOR(label, , setSizeFieldMapper)


forAllVectorNTypes(DEFINE_FIELD_MAPPER_OPERATOR, setSizeFieldMapper)

forAllTensorNTypes(DEFINE_FIELD_MAPPER_OPERATOR, setSizeFieldMapper)

forAllDiagTensorNTypes(DEFINE_FIELD_MAPPER_OPERATOR, setSizeFieldMapper)

forAllSphericalTensorNTypes(DEFINE_FIELD_MAPPER_OPERATOR, setSizeFieldMapper)


// ************************************************************************* //
