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
    (c) 2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::correctedCellVolumeWeightMethod

Description
    Cell-volume-weighted mesh-to-mesh interpolation class

    Volume conservative with calculated offset vectors

SourceFiles
    correctedCellVolumeWeightMethod.C

\*---------------------------------------------------------------------------*/

#ifndef correctedCellVolumeWeightMethod_H
#define correctedCellVolumeWeightMethod_H

#include "meshToMesh/calcMethod/cellVolumeWeight/cellVolumeWeightMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class correctedCellVolumeWeightMethod Declaration
\*---------------------------------------------------------------------------*/

class correctedCellVolumeWeightMethod
:
    public cellVolumeWeightMethod
{
protected:

    // Protected Member Functions

        //- Calculate the mesh-to-mesh addressing and weights
        void calculateAddressing
        (
            labelListList& srcToTgtCellAddr,
            scalarListList& srcToTgtCellWght,
            pointListList& srcToTgtCellVec,
            labelListList& tgtToSrcCellAddr,
            scalarListList& tgtToSrcCellWght,
            pointListList& tgtToSrcCellVec,
            const label srcSeedI,
            const label tgtSeedI,
            const labelList& srcCellIDs,
            boolList& mapFlag,
            label& startSeedI
        );


public:

    //- Run-time type information
    TypeName("correctedCellVolumeWeight");


    // Constructors

        //- Construct from source and target meshes
        correctedCellVolumeWeightMethod
        (
            const polyMesh& src,
            const polyMesh& tgt,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        correctedCellVolumeWeightMethod
        (
            const correctedCellVolumeWeightMethod&
        ) = delete;


    //- Destructor
    virtual ~correctedCellVolumeWeightMethod();


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
        void operator=(const correctedCellVolumeWeightMethod&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
