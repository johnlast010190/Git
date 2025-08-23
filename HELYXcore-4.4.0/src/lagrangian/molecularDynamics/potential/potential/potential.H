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
    (c) 2011 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::potential

Description

SourceFiles
    potentialI.H
    potential.C

\*---------------------------------------------------------------------------*/

#ifndef potential_H
#define potential_H

#include "meshes/polyMesh/polyMesh.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "db/Time/Time.H"
#include "pairPotential/pairPotentialList/pairPotentialList.H"
#include "electrostaticPotential/electrostaticPotential.H"
#include "tetherPotential/tetherPotentialList/tetherPotentialList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class potential Declaration
\*---------------------------------------------------------------------------*/

class potential
{
    // Private data

        const polyMesh& mesh_;

        List<word> idList_;

        List<word> siteIdList_;

        label nPairPotIds_;

        scalar potentialEnergyLimit_;

        labelList removalOrder_;

        pairPotentialList pairPotentials_;

        tetherPotentialList tetherPotentials_;

        vector gravity_;


    // Private Member Functions

        void setSiteIdList(const dictionary& moleculePropertiesDict);

        void readPotentialDict();

        void readMdInitialiseDict
        (
            const IOdictionary& mdInitialiseDict,
            IOdictionary& idListDict
        );


public:

    // Constructors

        //- Construct from mesh reference
        potential(const polyMesh& mesh);

        //- Construct from mdInitialiseDict
        potential
        (
            const polyMesh& mesh,
            const IOdictionary& mdInitialiseDict,
            IOdictionary& idListDict
        );

        //- Disallow default bitwise copy construction
        potential(const potential&) = delete;


    //- Destructor
    ~potential();


    // Member Functions

        // Access

            inline label nIds() const;

            inline const List<word>& idList() const;

            inline const List<word>& siteIdList() const;

            inline scalar potentialEnergyLimit() const;

            inline label nPairPotentials() const;

            inline const labelList& removalOrder() const;

            inline const pairPotentialList& pairPotentials() const;

            inline const tetherPotentialList& tetherPotentials() const;

            inline const vector& gravity() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const potential&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "potential/potentialI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
