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
    (c) 2011-2014 OpenFOAM Foundation
    (c) 2016-2024 Engys Ltd.

Class
    Foam::fvMeshTopoChangers::multiRefiner

Description
    A fvMesh with built-in multi-criterion refinement.

    Determines which cells to refine/unrefine and does all in update().
    Uses functionality of fvMeshTopoChangers::refiner.

SourceFiles
    fvMeshTopoChangersMultiRefiner.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshTopoChangersMultiRefiner_H
#define fvMeshTopoChangersMultiRefiner_H

#include "refiner/fvMeshTopoChangersRefiner.H"

namespace Foam
{
namespace fvMeshTopoChangers
{

/*---------------------------------------------------------------------------*\
               Class fvMeshTopoChangers::multiRefiner Declaration
\*---------------------------------------------------------------------------*/

class multiRefiner
:
    public refiner
{
    // Private Member Data

        //- Get access to refinement level and refinement field
        autoPtr<volScalarField> internalRefinementFieldPtr_;
        autoPtr<volScalarField> targetLevelPtr_;
        autoPtr<volScalarField> isLevelPtr_;

        //- Refinement based on field min/max
        HashTable< dictionary > fields_;

        //- Refinement based on mag(gradFields) min/max
        HashTable< List<scalar>> gradFields_;

        //- Refinement based on curl(field) min/max
        HashTable< List<scalar>> curlFields_;

        //- Refinement based on interface (multiphase solver)
        HashTable< dictionary > interface_;

        //- Refinement based on prescribed region
        PtrList<entry> refinedRegions_;

        //- Switch for standard refinement vs multi-criterion
        Switch enableRefinementControl_;

        //- Switch for multi-criterion refinement based on fixed target field
        Switch fixedTargetField_;

        //- Switch for target level smoothing
        Switch targetLevelSmoothing_;

        //- Switch to (dis-)allow target level growth
        Switch allowOnlyLevelGrowth_;

        //- Non-dimensionalized diffusion coefficient for elliptical relaxation
        scalar diffCoeff_;

        //- Smooth target refinement field only in vicinity of refinement jumps
        //  value. Controls how many neighbour cells are affected.
        label nLayerSmooth_;

        label nBufferLayers_;
        label refineInterval_;
        label maxRefinement_;


    // Private Member Functions

        //- (Re-)Calculate the refinement field
        void updateRefinementField();

        //- Read dynamicMeshDict
        void readDict();

        //- Get point min of connected cell
        virtual scalarField minCellField(const volScalarField&) const override;

        //- Overwrite protected member function of fvMeshTopoChangers::refiner
        virtual void selectRefineCandidates
        (
            const scalar lowerRefineLevel,
            const scalar upperRefineLevel,
            const scalarField& vFld,
            PackedBoolList& candidateCell
        ) const override;

        //- Smoothing algorithm for target field
        void smoothTargetLevel
        (
           const volScalarField& targetLevelField,
           labelList& targetLevel
        );


public:

    //- Runtime type information
    TypeName("multiRefiner");


    // Constructors

        //- Construct from fvMesh and dictionary
        multiRefiner(fvMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        multiRefiner(const multiRefiner&) = delete;


    //- Destructor
    virtual ~multiRefiner();


    // Member Functions

        //- Update the mesh for topology changes
        virtual bool update() override;

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&) override;

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&) override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiRefiner&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshTopoChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
