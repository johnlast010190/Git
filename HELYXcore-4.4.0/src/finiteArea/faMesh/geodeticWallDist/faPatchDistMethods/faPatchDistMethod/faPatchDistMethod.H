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
    (c) 2016-2024 Engys Ltd.
    (c) 2015 OpenFOAM Foundation

Class
    Foam::faPatchDistMethod

Description
    Specialisation of faPatchDist for geodetic wall distance calculation

SourceFiles
    faPatchDistMethod.C

\*---------------------------------------------------------------------------*/
#ifndef faPatchDistMethod_H
#define faPatchDistMethod_H

#include "db/dictionary/dictionary.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "fields/areaFields/areaFieldsFwd.H"
#include "meshes/polyMesh/polyTopoChangeMap/polyTopoChangeMap.H"
#include "fields/faPatchFields/basic/fixedValue/fixedValueFaPatchFields.H"
#include "fields/faPatchFields/basic/zeroGradient/zeroGradientFaPatchFields.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class faMesh;

/*---------------------------------------------------------------------------*\
                          Class faPatchDistMethod Declaration
\*---------------------------------------------------------------------------*/

class faPatchDistMethod
{
protected:

    // Protected Member Data

        //- Reference to the mesh
        const faMesh& mesh_;

        //- Set of patch IDs
        const labelHashSet patchIDs_;


public:

    //- Runtime type information
    TypeName("faPatchDistMethod");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            faPatchDistMethod,
            dictionary,
            (
                const dictionary& dict,
                const faMesh& mesh,
                const labelHashSet& patchIDs
            ),
            (dict, mesh, patchIDs)
        );


    // Constructors

        //- Construct from mesh and patch ID set
        faPatchDistMethod
        (
            const faMesh& mesh,
            const labelHashSet& patchIDs
        );

        //- Disallow default bitwise copy construction
        faPatchDistMethod(const faPatchDistMethod&) = delete;


    // Selectors

        static autoPtr<faPatchDistMethod> New
        (
            const dictionary& dict,
            const faMesh& mesh,
            const labelHashSet& patchIDs
        );


    //- Destructor
    virtual ~faPatchDistMethod();


    // Static Functions

        //- Return the patch types for y and n
        //  These are fixedValue for the set provided otherwise zero-gradient
        template<class Type>
        static wordList patchTypes
        (
            const faMesh& mesh,
            const labelHashSet& patchIDs
        );


    // Member Functions

        //- Return the patchIDs
        const labelHashSet& patchIDs() const
        {
            return patchIDs_;
        }

        //- Update cached geometry when the mesh moves
        virtual bool movePoints()
        {
            return true;
        }

        //- Update cached topology and geometry when the mesh changes
        virtual void topoChange(const polyTopoChangeMap&)
        {}

        //- Correct the given distance-to-patch field
        virtual bool correct(areaScalarField& y) = 0;

        //- Correct the given distance-to-patch and normal-to-patch fields
        virtual bool correct(areaScalarField& y, areaVectorField& n) = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const faPatchDistMethod&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "faMesh/geodeticWallDist/faPatchDistMethods/faPatchDistMethod/faPatchDistMethodTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
