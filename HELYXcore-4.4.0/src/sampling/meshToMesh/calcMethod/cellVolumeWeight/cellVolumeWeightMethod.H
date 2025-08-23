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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::cellVolumeWeightMethod

Description
    Cell-volume-weighted mesh-to-mesh interpolation class

    Volume conservative.

SourceFiles
    cellVolumeWeightMethod.C

\*---------------------------------------------------------------------------*/

#ifndef cellVolumeWeightMethod_H
#define cellVolumeWeightMethod_H

#include "meshToMesh/calcMethod/meshToMeshMethod/meshToMeshMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class cellVolumeWeightMethod Declaration
\*---------------------------------------------------------------------------*/

class cellVolumeWeightMethod
:
    public meshToMeshMethod
{
protected:

    // Protected Member Functions

        //- Find indices of overlapping cells in src and tgt meshes - returns
        //  true if found a matching pair
        bool findInitialSeeds
        (
            const labelList& srcCellIDs,
            const boolList& mapFlag,
            const label startSeedI,
            label& srcSeedI,
            label& tgtSeedI
        ) const;

        //- Calculate the mesh-to-mesh addressing and weights
        void calculateAddressing
        (
            labelListList& srcToTgtCellAddr,
            scalarListList& srcToTgtCellWght,
            labelListList& tgtToSrcCellAddr,
            scalarListList& tgtToSrcCellWght,
            const label srcSeedI,
            const label tgtSeedI,
            const labelList& srcCellIDs,
            boolList& mapFlag,
            label& startSeedI
        );

        //- Set the next cells in the advancing front algorithm
        void setNextCells
        (
            label& startSeedI,
            label& srcCelli,
            label& tgtCelli,
            const labelList& srcCellIDs,
            const boolList& mapFlag,
            const DynamicList<label>& visitedCells,
            labelList& seedCells
        ) const;


public:

    //- Run-time type information
    TypeName("cellVolumeWeight");


    // Constructors

        //- Construct from source and target meshes
        cellVolumeWeightMethod
        (
            const polyMesh& src,
            const polyMesh& tgt,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        cellVolumeWeightMethod(const cellVolumeWeightMethod&) = delete;


    //- Destructor
    virtual ~cellVolumeWeightMethod();


    // Member Functions

        //- Calculate addressing and weights and optionally offset vectors
        virtual void calculate
        (
            labelListList& srcToTgtAddr,
            scalarListList& srcToTgtWght,
            pointListList& srcToTgtVec,
            labelListList& tgtToSrcAddr,
            scalarListList& tgtToSrcWght,
            pointListList& tgtToSrcVec
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cellVolumeWeightMethod&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
