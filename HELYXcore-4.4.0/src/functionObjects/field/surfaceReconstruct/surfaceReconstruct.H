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
    (c) 2016 OpenCFD Ltd.
    (c) 2019 Engys Ltd.

Class
    Foam::functionObjects::surfaceReconstruct

Group
    grpFieldFunctionObjects

Description
    Linearly reconstruct surface fields to generate volume fields.

    Fields are stored
    - every time step the field is updated with new values
    - at output it writes the fields

    This functionObject can either be used
    - to calculate a new field as a  post-processing step or
    - since the fields are registered, used in another functionObject

Usage
    Example of function object specification:
    \verbatim
    surfaceReconstruct1
    {
        type        surfaceReconstruct;
        libs        ("libfieldFunctionObjects.so");
        ...
        fields      ((p pNear)(U UNear));
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property | Description               | Required    | Default value
        type     | type name: nearWallFields | yes         |
        fields   | list of fields with correspoding output field names | yes |
        log      | Log to standard output    | no          | yes
    \endtable


See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::timeControl

SourceFiles
    surfaceReconstruct.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_surfaceReconstruct_H
#define functionObjects_surfaceReconstruct_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "primitives/Tuple2/Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyTopoChangeMap;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                  Class surfaceReconstruct Declaration
\*---------------------------------------------------------------------------*/

class surfaceReconstruct
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- Fields to process
        List<Tuple2<word, word>> fieldSet_;

        //- Switch to send output to Info as well as to file
        Switch log_;


    // Protected Member Functions

        template<class Type>
        void reconstructFields();


public:

    //- Runtime type information
    TypeName("surfaceReconstruct");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        surfaceReconstruct
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        surfaceReconstruct(const surfaceReconstruct&) = delete;


    //- Destructor
    virtual ~surfaceReconstruct();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Calculate the reconstruct fields
        virtual bool execute();

        //- Write the reconstruct fields
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const surfaceReconstruct&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfaceReconstruct/surfaceReconstructTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
