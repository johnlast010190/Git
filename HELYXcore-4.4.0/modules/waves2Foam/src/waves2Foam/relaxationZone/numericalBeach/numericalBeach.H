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
    (c) held by original author
    (c) 2024 Engys Ltd.

Class
    Foam::numericalBeaches::numericalBeach

Description
    An abstract base class for different types of numerical beaches

SourceFiles
    numericalBeach.C
    newNumericalBeach.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef numericalBeach_H
#define numericalBeach_H

#include "fields/volFields/volFields.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace numericalBeaches
{

/*---------------------------------------------------------------------------*\
                      Class numericalBeach Declaration
\*---------------------------------------------------------------------------*/

class numericalBeach
:
    public IOdictionary
{
protected:

    // Protected data

        const fvMesh& mesh_;

        const dictionary coeffDict_;

        word beachType_;


public:

    //- Runtime type information
    TypeName("numericalBeach");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            numericalBeach,
            dictionary,
            (
                const word& subDictName, const fvMesh& mesh_
            ),
            (subDictName, mesh_)
        );


    // Constructors

        //- Construct from components
        numericalBeach
        (
            const word& type,
            const fvMesh& mesh_
        );

        //- Disallow default bitwise copy construction
        numericalBeach(const numericalBeach&) = delete;


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<numericalBeach> New
        (
            const word& subDictName,
            const fvMesh& mesh_
        );


    //- Destructor
    virtual ~numericalBeach();


    // Member Functions

        virtual void correct
        (
            const labelList&,
            const scalarField&,
            volScalarField&
        ) = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const numericalBeach&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationShapes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
