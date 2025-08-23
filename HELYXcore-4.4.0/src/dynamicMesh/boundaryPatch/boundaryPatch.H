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
    (c) 2017-2022 Engys Ltd.

Class
    Foam::boundPatch

Description
    Like polyPatch but without reference to mesh. patchIdentifier::index
    is not used. Used in boundaryMesh to hold data on patches.

SourceFiles
    boundaryPatch.C

\*---------------------------------------------------------------------------*/

#ifndef boundPatch_H
#define boundPatch_H

#include "meshes/Identifiers/patch/patchIdentifier.H"
#include "memory/autoPtr/autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class boundPatch;

Ostream& operator<<(Ostream&, const boundPatch&);


/*---------------------------------------------------------------------------*\
                           Class boundPatch Declaration
\*---------------------------------------------------------------------------*/

class boundPatch
:
    public patchIdentifier
{
    // Private data

        label size_;
        label start_;

public:

    // Constructors

        //- Construct from components
        boundPatch
        (
            const word& name,
            const label index,
            const label size,
            const label start,
            const word& physicalType = word::null
        );

        //- Construct from dictionary
        boundPatch
        (
            const word& name,
            const dictionary& dict,
            const label index
        );

        //- Construct as copy, resetting the index
        boundPatch(const boundPatch&, const label index);

        //- Clone
        autoPtr<boundPatch> clone() const;


    //- Destructor
    ~boundPatch();


    // Member Functions

        label size() const
        {
            return size_;
        }

        label& size()
        {
            return size_;
        }

        label start() const
        {
            return start_;
        }

        label& start()
        {
            return start_;
        }


        //- Write dictionary
        virtual void write(Ostream&) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const boundPatch&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
