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
    (c) 2016 OpenCFD Ltd.
    (c) 2013-2018 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::mapNearestAMI

Description
    Nearest-mapping Arbitrary Mesh Interface (AMI) method

SourceFiles
    mapNearestAMI.C

\*---------------------------------------------------------------------------*/

#ifndef mapNearestAMI_H
#define mapNearestAMI_H

#include "AMIInterpolation/AMIInterpolation/AMIMethod/AMIMethod/AMIMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class mapNearestAMI Declaration
\*---------------------------------------------------------------------------*/

class mapNearestAMI
:
    public AMIMethod
{
private:

    // Private Member Functions

        // Marching front

            //- Find nearest target face for source face srcFacei
            void findNearestFace
            (
                const primitivePatch& srcPatch,
                const primitivePatch& tgtPatch,
                const label& srcFacei,
                label& tgtFacei
            ) const;

            //- Determine next source-target face pair
            void setNextNearestFaces
            (
                boolList& mapFlag,
                label& startSeedI,
                label& srcFacei,
                label& tgtFacei
            ) const;

            //- Find mapped source face
            label findMappedSrcFace
            (
                const label tgtFacei,
                const List<DynamicList<label>>& tgtToSrc
            ) const;


        // Evaluation

            //- Area of intersection between source and target faces
            scalar interArea
            (
                const label srcFacei,
                const label tgtFacei
            ) const;


public:

    //- Runtime type information
    TypeName("mapNearestAMI");


    // Constructors

        //- Construct from components
        mapNearestAMI
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
        mapNearestAMI(const mapNearestAMI&) = delete;


    //- Destructor
    virtual ~mapNearestAMI();


    // Member Functions

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

            //- Normalise the weight. Can optionally subset addressing
            //  (e.g. for mapNearest)
            virtual void normaliseWeights
            (
                const bool verbose,
                AMIInterpolation& inter
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const mapNearestAMI&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
