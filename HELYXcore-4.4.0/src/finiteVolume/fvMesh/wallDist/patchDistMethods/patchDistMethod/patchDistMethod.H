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
    (c) 2017 OpenCFD Ltd.
    (c) 2015-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::patchDistMethod

Description
    Specialisation of patchDist for wall distance calculation

SourceFiles
    patchDistMethod.C

\*---------------------------------------------------------------------------*/

#ifndef patchDistMethod_H
#define patchDistMethod_H

#include "db/dictionary/dictionary.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "fields/volFields/volFieldsFwd.H"
#include "meshes/polyMesh/polyTopoChangeMap/polyTopoChangeMap.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                          Class patchDistMethod Declaration
\*---------------------------------------------------------------------------*/

class patchDistMethod
{
protected:

    // Protected Member Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Set of patch IDs
        const labelHashSet patchIDs_;


public:

    //- Runtime type information
    TypeName("patchDistMethod");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            patchDistMethod,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh,
                const labelHashSet& patchIDs
            ),
            (dict, mesh, patchIDs)
        );


    // Constructors

        //- Construct from mesh and patch ID set
        patchDistMethod
        (
            const fvMesh& mesh,
            const labelHashSet& patchIDs
        );

        //- Disallow default bitwise copy construction
        patchDistMethod(const patchDistMethod&) = delete;


    // Selectors

        static autoPtr<patchDistMethod> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const labelHashSet& patchIDs
        );

        static autoPtr<patchDistMethod> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const labelHashSet& patchIDs,
            const word& defaultPatchDistMethod
        );


    //- Destructor
    virtual ~patchDistMethod();


    // Static Functions

        //- Return the patch types for y and n
        //  These are fixedValue for the set provided otherwise zero-gradient
        template<class Type>
        static wordList patchTypes
        (
            const fvMesh& mesh,
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

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&)
        {}

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&)
        {}

        //- Correct the given distance-to-patch field
        virtual bool correct(volScalarField& y) = 0;

        //- Correct the given distance-to-patch and normal-to-patch fields
        virtual bool correct(volScalarField& y, volVectorField& n) = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const patchDistMethod&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMesh/wallDist/patchDistMethods/patchDistMethod/patchDistMethodTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
