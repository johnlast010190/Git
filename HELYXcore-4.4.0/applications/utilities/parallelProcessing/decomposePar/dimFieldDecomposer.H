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
    (c) 2022-2024 Engys Ltd.

Class
    Foam::dimFieldDecomposer

Description
    Dimensioned field decomposer.

SourceFiles
    dimFieldDecomposer.C
    dimFieldDecomposerDecomposeFields.C

\*---------------------------------------------------------------------------*/

#ifndef dimFieldDecomposer_H
#define dimFieldDecomposer_H

#include "fvMesh/fvMesh.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOobjectList;

/*---------------------------------------------------------------------------*\
                      Class dimFieldDecomposer Declaration
\*---------------------------------------------------------------------------*/

class dimFieldDecomposer
{
    // Private Data

        //- Reference to processor mesh
        const fvMesh& procMesh_;

        //- Reference to cell addressing
        const labelList& cellAddressing_;


public:

    // Constructors

        //- Construct from components
        dimFieldDecomposer
        (
            const fvMesh& completeMesh,
            const fvMesh& procMesh,
            const labelList& faceAddressing,
            const labelList& cellAddressing
        );

        //- Disallow default bitwise copy construction
        dimFieldDecomposer(const dimFieldDecomposer&) = delete;


    //- Destructor
    ~dimFieldDecomposer();


    // Member Functions

        //- Decompose field
        template<class Type>
        tmp<DimensionedField<Type, volMesh>> decomposeField
        (
            const DimensionedField<Type, volMesh>& field
        ) const;

        //- Decompose llist of fields
        template<class GeoField>
        void decomposeFields(const PtrList<GeoField>& fields) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dimFieldDecomposer&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dimFieldDecomposerDecomposeFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
