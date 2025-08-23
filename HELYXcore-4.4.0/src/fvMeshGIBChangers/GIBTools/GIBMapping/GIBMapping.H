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

Class
    Foam::GIBMapping

Description
    Toolbox for GIB mapping.

SourceFiles
    GIBMapping.C

\*---------------------------------------------------------------------------*/

#ifndef GIBMapping_H
#define GIBMapping_H

#include "base/fvMeshGIBChangersBase.H"
#include "AMIInterpolation/AMIInterpolation/AMIInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMeshGIBChangersBase;

/*---------------------------------------------------------------------------*\
                            Class GIBMapping Declaration
\*---------------------------------------------------------------------------*/

class GIBMapping
{
    // Private Data

        const fvMeshGIBChangersBase& GIBChanger_;

        bool mapSlave_;

        //- AMI matching settings
        const scalar AMICosMatchAngle_;
        const scalar maxAMIWeightScale_;

        List<Tuple2< word, Field<scalar >>> gmfvScalars_;
        List<Tuple2< word, Field<vector>>> gmfvVectors_;
        List<Tuple2< word, Field<sphericalTensor>>> gmfvSpTensors_;
        List<Tuple2< word, Field<symmTensor>>> gmfvSymmTensors_;
        List<Tuple2< word, Field<tensor>>> gmfvTensors_;

        List<Tuple2< word, Field<scalar>>> gmfvsScalars_;
        List<Tuple2< word, Field<vector>>> gmfvsVectors_;
        List<Tuple2< word, Field<sphericalTensor>>> gmfvsSpTensors_;
        List<Tuple2< word, Field<symmTensor>>> gmfvsSymmTensors_;
        List<Tuple2< word, Field<tensor>>> gmfvsTensors_;

        List<Tuple2< word, Field<scalar>>> gsfvScalars_;
        List<Tuple2< word, Field<vector>>> gsfvVectors_;
        List<Tuple2< word, Field<sphericalTensor>>> gsfvSpTensors_;
        List<Tuple2< word, Field<symmTensor>>> gsfvSymmTensors_;
        List<Tuple2< word, Field<tensor>>> gsfvTensors_;

        List<Tuple2< word, Field<scalar>>> gsfvsScalars_;
        List<Tuple2< word, Field<vector>>> gsfvsVectors_;
        List<Tuple2< word, Field<sphericalTensor>>> gsfvsSpTensors_;
        List<Tuple2< word, Field<symmTensor>>> gsfvsSymmTensors_;
        List<Tuple2< word, Field<tensor>>> gsfvsTensors_;

        mutable AMIInterpolation* mAMIInterPtr_;
        mutable AMIInterpolation* sAMIInterPtr_;

        const primitivePatch& sourcePatch_;

        //- Local patch data
        mutable primitivePatch* sourcePatchm_;
        mutable primitivePatch* sourcePatchs_;

        Field<vector> sourcePointsm_;
        Field<vector> sourcePointss_;
        faceList sourceFacesm_;
        faceList sourceFacess_;

        mutable primitivePatch* targetPatchm_;
        mutable primitivePatch* targetPatchs_;
        Field<vector> targetPointsm_;
        Field<vector> targetPointss_;
        faceList targetFacesm_;
        faceList targetFacess_;

        //- Global patch data
        Field<vector> gppPoints_;
        faceList gppFaces_;


    // Private Member Functions

        void storeFields();

        void mapSerial();

        void combinePolyPatch();

        void mapMaster();

        void mapSlave();

        void mapNonOverlappingTargetFaces();

        labelList calcNearestPatchFaceMapping
        (
            const labelList&,
            const primitivePatch&,
            const primitivePatch&,
            const labelList&
        ) const;

        void makeMasterAMIInterpolation() const;

        void makeSlaveAMIInterpolation() const;

        //- Clear all demand-driven data
        virtual void clearOutGIBData();

        template<class Type, template<class> class PatchField, class GeoMesh>
        void StoreOldFieldsToPatch
        (
            List< Tuple2<word, Field<Type>>>&,
            const label&
        );

        template<class Type, template<class> class PatchField, class GeoMesh>
        void MapGIBField
        (
            List<Tuple2<word, Field<Type>>>&,
            const label&,
            const AMIInterpolation&
        );

        template<class Type, template<class> class PatchField, class GeoMesh>
        void Map1to1GIBField
        (
            const List<Tuple2<word, Field<Type>>>&,
            const List<Tuple2<word, Field<Type>>>&
        );

        template<class Type, template<class> class PatchField, class GeoMesh>
        void MapNonOverlapFaces
        (
            const List<Tuple2<word, Field<Type>>>&,
            const label,
            const labelList&,
            const labelList&
        );


public:

    //- Runtime type information
    TypeName("GIBMapping");


    // Constructors

        //- Construct from components
        GIBMapping
        (
            const fvMeshGIBChangersBase&,
            const primitivePatch&,
            bool mapSlave = true,
            // Default to the most permissive possible AMI settings in the GIB
            const scalar& AMICosMatchAngle = -1,
            const scalar& maxAMIWeightScale = GREAT
        );


    //- Destructor
    virtual ~GIBMapping()
    {
        clearOutGIBData();
    };


    // Member Functions

        virtual void mapBcs();

        virtual void map1To1();

        virtual void mapNonOverlapFaces
        (
            const labelList&,
            const labelList&,
            const label
        );

        virtual triSurface triS();

        const AMIInterpolation& masterAMIInter() const
        {
            if (!mAMIInterPtr_)
            {
                makeMasterAMIInterpolation();
            }
            return *mAMIInterPtr_;
        }

        const AMIInterpolation& slaveAMIInter() const
        {
            if (!sAMIInterPtr_)
            {
                makeSlaveAMIInterpolation();
            }
            return *sAMIInterPtr_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "GIBTools/GIBMapping/GIBMappingTemplates.C"
#endif

#endif

// ************************************************************************* //
