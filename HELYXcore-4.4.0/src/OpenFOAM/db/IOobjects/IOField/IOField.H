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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2018 OpenCFD Ltd.
    (c) 2022 Engys Ltd.

Class
    Foam::IOField

Description
    A primitive field of type \<T\> with automated input and output.

SourceFiles
    IOField.C

\*---------------------------------------------------------------------------*/

#ifndef IOField_H
#define IOField_H

#include "db/regIOobject/regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IOField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class IOField
:
    public regIOobject,
    public Field<Type>
{

public:

    TypeName("Field");


    // Constructors

        //- Default copy construct
        IOField(const IOField&) = default;

        //- Construct from IOobject
        IOField(const IOobject&);

        //- Construct from IOobject; does local processor require reading?
        IOField(const IOobject&, const bool read);

        //- Construct from IOobject and size (does not set values)
        IOField(const IOobject&, const label size);

        //- Construct from components
        IOField(const IOobject&, const Field<Type>&);

        //- Construct by transferring the Field contents
        IOField(const IOobject&, Field<Type>&&);

        //- Move construct by transferring the Field contents
        IOField(const IOobject&, const tmp<Field<Type>>&);

        //- Move constructor
        IOField(IOField<Type>&&);


    //- Destructor
    virtual ~IOField();


    // Member functions

        bool writeData(Ostream&) const;


    // Member operators

        void operator=(const IOField<Type>&);
        void operator=(IOField<Type>&&);

        void operator=(const Field<Type>&);
        void operator=(Field<Type>&&);
        void operator=(const tmp<Field<Type>>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "db/IOobjects/IOField/IOField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
