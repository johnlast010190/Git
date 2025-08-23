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
    (c) 2016-2020 OpenCFD Ltd.
    (c) 2013-2018 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::partialFaceAreaWeightAMI

Description
    Partial face area weighted Arbitrary Mesh Interface (AMI) method

SourceFiles
    partialFaceAreaWeightAMI.C

\*---------------------------------------------------------------------------*/

#ifndef partialFaceAreaWeightAMI_H
#define partialFaceAreaWeightAMI_H

#include "AMIInterpolation/AMIInterpolation/AMIMethod/faceAreaWeightAMI/faceAreaWeightAMI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class partialFaceAreaWeightAMI Declaration
\*---------------------------------------------------------------------------*/

class partialFaceAreaWeightAMI
:
    public faceAreaWeightAMI
{

private:

    // Private Member Functions

        // Marching front

            //- Set the source and target seed faces
            virtual void setNextFaces
            (
                label& startSeedI,
                label& srcFacei,
                label& tgtFacei,
                const boolList& mapFlag,
                labelList& seedFaces,
                const DynamicList<label>& visitedFaces,
                bool errorOnNotFound = true
            ) const;


public:

    //- Runtime type information
    TypeName("partialFaceAreaWeightAMI");


    // Constructors

        //- Construct from components
        partialFaceAreaWeightAMI
        (
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch,
            const scalarField& srcMagSf,
            const scalarField& tgtMagSf,
            const faceAreaIntersect::triangulationMode& triMode,
            const scalar& cosMatchAngle,
            const bool reverseTarget = false,
            const bool requireMatch = true
        );

        //- Disallow default bitwise copy construction
        partialFaceAreaWeightAMI(const partialFaceAreaWeightAMI&) = delete;


    //- Destructor
    virtual ~partialFaceAreaWeightAMI();


    // Member Functions

        // Access

            //- Flag to indicate that interpolation patches are conformal
            virtual bool conformal() const;


        // Manipulation

            //- Update addressing and weights
            virtual void calculate
            (
                labelListList& srcAddress,
                scalarListList& srcWeights,
                scalarList& symSrcWeights,
                pointListList& srcCentroids,
                labelListList& tgtAddress,
                scalarListList& tgtWeights,
                scalarList& symTgtWeights,
                label srcFacei = -1,
                label tgtFacei = -1
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const partialFaceAreaWeightAMI&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
