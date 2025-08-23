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
    (c) 2016 OpenCFD Ltd.

Class
    Foam::ensightSerialOutput

Description
    A collection of functions for writing ensight file content in serial.
    Can be considered transitional and may eventually merge with ensightOutput.

SourceFiles
    ensightSerialOutputTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef ensightSerialOutput_H
#define ensightSerialOutput_H

#include "ensight/part/ensightPart.H"
#include "ensight/part/ensightPartFaces.H"
#include "ensight/part/ensightPartCells.H"
#include "ensight/part/ensightParts.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ensightSerialOutput Declaration
\*---------------------------------------------------------------------------*/

class ensightSerialOutput
{
    // Private Methods

    //- Write field content (component-wise) for the given ensight element type
    template<template<typename> class FieldContainer, class Type>
    static void writeFieldContent
    (
        const word& key,
        const FieldContainer<Type>& fld,
        ensightFile& os
    );


    //- Disallow null constructor
    ensightSerialOutput() = delete;

public:

    // Public Methods

    //- Write field component-wise
    template<class Type>
    static bool writeField
    (
        const Field<Type>&,
        const ensightPartFaces& part,
        ensightFile& os,
        const bool nodeValues
    );


    //- Write volume field component-wise
    template<class Type>
    static bool writeField
    (
        const VolField<Type>&,
        const ensightPartCells& part,
        ensightFile& os
    );


    //- Write volume field component-wise
    template<class Type>
    static bool writeField
    (
        const VolField<Type>&,
        const ensightParts& list,
        ensightFile& os
    );


    //- Write volume field component-wise
    template<class Type>
    static inline bool writeField
    (
        const VolField<Type>& vf,
        const ensightPartCells& part,
        autoPtr<ensightFile>& output
    )
    {
        return writeField(vf, part, output.rawRef());
    }


    //- Write volume field component-wise
    template<class Type>
    static inline bool writeField
    (
        const VolField<Type>& vf,
        const ensightParts& list,
        autoPtr<ensightFile>& output
    )
    {
        return writeField(vf, list, output.rawRef());
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ensight/output/ensightSerialOutputTemplates.C"
#endif

#endif

// ************************************************************************* //
