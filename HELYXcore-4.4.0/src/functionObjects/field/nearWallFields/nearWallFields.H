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
    (c) 2015 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::nearWallFields

Group
    grpFieldFunctionObjects

Description
    Samples near-patch volume fields.

    Fields are stored
    - every time step the field is updated with new values
    - at output it writes the fields

    This functionObject can either be used
    - to calculate a new field as a  post-processing step or
    - since the fields are registered, used in another functionObject

Usage
    Example of function object specification:
    \verbatim
    nearWallFields1
    {
        type        nearWallFields;
        libs        ("libfieldFunctionObjects.so");

        writeControl writeTime;

        fields
        (
            (p pNear)
            (U UNear)
        );

        patches     (movingWall);

        distance    0.13;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property | Description               | Required    | Default value
        type     | type name: nearWallFields | yes         |
        fields   | list of fields with corresponding output field names | yes |
        patches  | list of patches to sample | yes         |
        distance | distance from patch to sample | yes     |
        log      | Log to standard output    | no          | yes
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    nearWallFields.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_nearWallFields_H
#define functionObjects_nearWallFields_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "fields/volFields/volFields.H"
#include "primitives/Tuple2/Tuple2.H"
#include "interpolation/interpolation/interpolationCellPoint/interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class nearWallFields Declaration
\*---------------------------------------------------------------------------*/

class nearWallFields
:
    public fvMeshFunctionObject
{
protected:

    // Protected member data

        // Read from dictionary

            //- Fields to process
            List<Tuple2<word, word>> fieldSet_;

            //- Switch to send output to Info as well as to file
            Switch log_;

            //- Patches to sample
            labelHashSet patchSet_;

            //- Distance away from wall
            scalar distance_;

            //- From original field to sampled result
            HashTable<word> fieldMap_;

            //- From resulting back to original field
            HashTable<word> reverseFieldMap_;


        // Calculated addressing

            //- From cell to seed patch faces
            labelListList cellToWalls_;

            //- From cell to tracked end point
            List<List<point>> cellToSamples_;

            //- Map from cell based data back to patch based data
            autoPtr<distributionMap> getPatchDataMapPtr_;


        // Locally constructed fields

            PtrList<volScalarField> vsf_;
            PtrList<volVectorField> vvf_;
            PtrList<volSphericalTensorField> vSpheretf_;
            PtrList<volSymmTensorField> vSymmtf_;
            PtrList<volTensorField> vtf_;


    // Protected Member Functions

        //- Calculate addressing from cells back to patch faces
        void calcAddressing();

        template<class Type>
        void createFields(PtrList<VolField<Type>>&) const;

        //- Override boundary fields with sampled values
        template<class Type>
        void sampleBoundaryField
        (
            const interpolationCellPoint<Type>& interpolator,
            VolField<Type>& fld
        ) const;

        template<class Type>
        void sampleFields(PtrList<VolField<Type>>&) const;


public:

    //- Runtime type information
    TypeName("nearWallFields");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        nearWallFields
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        nearWallFields(const nearWallFields&) = delete;


    //- Destructor
    virtual ~nearWallFields();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Calculate the near-wall fields
        virtual bool execute();

        //- Write the near-wall fields
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const nearWallFields&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nearWallFields/nearWallFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
