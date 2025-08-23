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
    (c) 2016 OpenFOAM Foundation

Class
    Foam::functionObjects::add

Group
    grpFieldFunctionObjects

Description
    Add a list of fields.

    The operation can be applied to any volume or surface fields generating a
    volume or surface scalar field.

    Example of function object specification:
    \verbatim
    Ttot
    {
        type            add;
        libs            ("libfieldFunctionObjects.so");
        fields          (T Tdelta);
        result          Ttot;
        executeControl  writeTime;
        writeControl    writeTime;
    }
    \endverbatim

See also
    Foam::functionObjects::fieldsExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    add.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_add_H
#define functionObjects_add_H

#include "fieldsExpression/fieldsExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class add Declaration
\*---------------------------------------------------------------------------*/

class add
:
    public fieldsExpression
{
    // Private Member Functions

        //- Add the list of fields of the specified type
        //  and return the result
        template<class GeoFieldType>
        tmp<GeoFieldType> calcFieldType() const;

        //- Add the list of fields and return true if successful
        virtual bool calc();


public:

    friend class fieldsExpression;


    //- Runtime type information
    TypeName("add");


    // Constructors

        //- Construct from Time and dictionary
        add
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~add();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "add/addTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
