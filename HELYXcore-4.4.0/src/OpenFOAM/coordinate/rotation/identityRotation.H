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
    (c) 2018 OpenCFD Ltd.

Class
    Foam::coordinateRotations::identity

Description
    An identity coordinateRotation.

    \verbatim
    coordinateRotation
    {
        type    none;
    }
    \endverbatim

    \heading Dictionary entries
    \table
        Property    | Description                           | Required | Default
        type        | Type name: none                       | yes      |
    \endtable

SourceFiles
    identityRotation.C

\*---------------------------------------------------------------------------*/

#ifndef coordinateRotations_identity_H
#define coordinateRotations_identity_H

#include "coordinate/rotation/coordinateRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coordinateRotations
{

/*---------------------------------------------------------------------------*\
                Class coordinateRotations::identity Declaration
\*---------------------------------------------------------------------------*/

class identity
:
    public coordinateRotation
{
public:

    //- Runtime type information
    TypeNameNoDebug("none");


    // Constructors

        //- Construct null
        identity();

        //- Copy construct
        identity(const identity& unused);

        //- Construct from dictionary
        explicit identity(const dictionary& unused);

        //- Return clone
        autoPtr<coordinateRotation> clone() const
        {
            return
                autoPtr<coordinateRotation>::NewFrom
                <coordinateRotations::identity>(*this);
        }


    //- Destructor
    virtual ~identity() = default;


    // Member Functions

        //- Reset specification (no-op)
        virtual void clear();

        //- Return an identity rotation tensor
        virtual tensor R() const;

        //- Update rotation
        virtual void update
        (
            const vector& e1,
            const vector& e2
        ){};

        //- Write information
        virtual void write(Ostream& os) const;

        //- Write dictionary entry
        virtual void writeEntry(const word& keyword, Ostream& os) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coordinateRotations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
