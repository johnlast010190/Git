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

SourceFiles
    tableReaders.C

\*---------------------------------------------------------------------------*/

#ifndef tableReaders_H
#define tableReaders_H

#include "interpolations/interpolationTable/tableReaders/tableReader.H"
#include "fields/Fields/fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Only used internally
#define makeTypeTableReadersTypeName(typeTableReader, dataType)                \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(typeTableReader<dataType>, 0)

// Sometimes used externally
#define makeTableReadersTypeName(typeTableReader)                              \
                                                                               \
    makeTypeTableReadersTypeName(typeTableReader, scalar);                     \
    makeTypeTableReadersTypeName(typeTableReader, vector);                     \
    makeTypeTableReadersTypeName(typeTableReader, sphericalTensor);            \
    makeTypeTableReadersTypeName(typeTableReader, symmTensor);                 \
    makeTypeTableReadersTypeName(typeTableReader, tensor)

// Define type info for single dataType template instantiation (eg, vector)
#define makeTableReaderType(typeTableReader, dataType)                         \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(typeTableReader<dataType>, 0);         \
    addTemplatedToRunTimeSelectionTable                                        \
    (                                                                          \
        tableReader, typeTableReader, dataType, dictionary                     \
    )


// Define type info for scalar, vector etc. instantiations
#define makeTableReaders(typeTableReader)                                      \
                                                                               \
    makeTableReaderType(typeTableReader, scalar);                              \
    makeTableReaderType(typeTableReader, vector);                              \
    makeTableReaderType(typeTableReader, sphericalTensor);                     \
    makeTableReaderType(typeTableReader, symmTensor);                          \
    makeTableReaderType(typeTableReader, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
