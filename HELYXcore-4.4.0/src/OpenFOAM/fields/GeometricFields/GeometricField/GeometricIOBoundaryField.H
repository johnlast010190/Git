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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::GeometricIOBoundaryField

Description
    GeometricBoundaryField class with the ability to store on the object
    registry. For storing standalone GeometricBoundaryFields.

SourceFiles
    GeometricIOBoundaryField.C

\*---------------------------------------------------------------------------*/

#ifndef GeometricIOBoundaryField_H
#define GeometricIOBoundaryField_H

#include "dimensionedTypes/dimensionedTypes.H"
#include "GeometricBoundaryField.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"
#include "fields/FieldFields/FieldField/FieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;

/*---------------------------------------------------------------------------*\
                    Class GeometricIOBoundaryField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricIOBoundaryField
:
    public regIOobject,
    public GeometricBoundaryField<Type, PatchField, GeoMesh>
{
public:

    // Public typedefs

        //- Parent class
        typedef GeometricBoundaryField<Type, PatchField, GeoMesh>
            GeoBoundaryField;

        //- Inherited typedef from parent
        typedef typename GeoBoundaryField::BoundaryMesh BoundaryMesh;

        //- Inherited typedef from parent
        typedef typename GeoBoundaryField::Internal Internal;

    // Constructors

        //- Construct from a BoundaryMesh
        GeometricIOBoundaryField(const IOobject& io, const BoundaryMesh&);

        //- Construct from a BoundaryMesh, reference to the internal field
        //  and a patch field type
        GeometricIOBoundaryField
        (
            const IOobject& io,
            const BoundaryMesh&,
            const Internal&,
            const word&,
            const wordList& actualPatchTypes = wordList()
        );

        //- Construct from a BoundaryMesh, reference to the internal field,
        //  a wordList of patch field types and optionally
        //  the actual patch types (to override constraint patches)
        GeometricIOBoundaryField
        (
            const IOobject& io,
            const BoundaryMesh&,
            const Internal&,
            const wordList& wantedPatchTypes,
            const wordList& actualPatchTypes = wordList()
        );

        //- Construct from a BoundaryMesh, reference to the internal field
        //  and a PtrList<PatchField<Type>>
        GeometricIOBoundaryField
        (
            const IOobject& io,
            const BoundaryMesh&,
            const Internal&,
            const PtrList<PatchField<Type>>&
        );

        //- Construct as copy setting the reference to the internal field
        GeometricIOBoundaryField
        (
            const IOobject& io,
            const Internal&,
            const GeoBoundaryField&
        );

        //- Construct as copy, or move, setting the reference to the internal field
        GeometricIOBoundaryField
        (
            const IOobject& io,
            const Internal&,
            GeoBoundaryField&,
            bool reuse
        );

        //- Construct as copy
        //  Dangerous because Field may be set to a field which gets deleted
        //  Need new type of BoundaryField, one which is part of a geometric
        //  field for which snGrad etc. may be called and a free standing
        //  BoundaryField for which such operations are unavailable.
        GeometricIOBoundaryField(const GeometricIOBoundaryField&);

        //- Construct from dictionary
        GeometricIOBoundaryField
        (
            const IOobject& io,
            const BoundaryMesh&,
            const Internal&,
            const dictionary&
        );

        // Public member functions

            // Write

            virtual bool writeData(Ostream& os) const
            {
                os << *this;
                return os.good();
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/GeometricFields/GeometricField/GeometricIOBoundaryField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
