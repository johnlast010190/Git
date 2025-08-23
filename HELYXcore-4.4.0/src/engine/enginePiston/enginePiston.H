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
    Foam::enginePiston

Description
    Foam::enginePiston

SourceFiles
    enginePiston.C

\*---------------------------------------------------------------------------*/

#ifndef enginePiston_H
#define enginePiston_H

#include "meshes/Identifiers/patch/polyPatchID.H"
#include "coordinate/systems/coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class engineTime;

/*---------------------------------------------------------------------------*\
                           Class enginePiston Declaration
\*---------------------------------------------------------------------------*/

class enginePiston
{
    // Private data

        //- Reference to engine mesh
        const polyMesh& mesh_;

        //- Reference to engine database
        const engineTime& engineDB_;

        //- Piston patch
        polyPatchID patchID_;

        //- Coordinate system
        autoPtr<coordinateSystem> csPtr_;


        // Piston layering data

            //- Min layer thickness
            const scalar minLayer_;

            //- Max layer thickness
            const scalar maxLayer_;


public:

    // Constructors

        //- Construct from components
        enginePiston
        (
            const polyMesh& mesh,
            const word& pistonPatchName,
            const autoPtr<coordinateSystem>& pistonCS,
            const scalar minLayer,
            const scalar maxLayer

        );

        //- Construct from dictionary
        enginePiston
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        enginePiston(const enginePiston&) = delete;


    // Member Functions

        //- Return coordinate system
        const coordinateSystem& cs() const
        {
            return csPtr_();
        }

        //- Return ID of piston patch
        const polyPatchID& patchID() const
        {
            return patchID_;
        }

        // Piston layering thickness

            scalar minLayer() const
            {
                return minLayer_;
            }

            scalar maxLayer() const
            {
                return maxLayer_;
            }


        //- Write dictionary
        void writeDict(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const enginePiston&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
