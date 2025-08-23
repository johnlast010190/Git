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
    (c) 2015 OpenCFD Ltd.
    (c) 2018-2020 OpenFOAM Foundation
    (c) 2021-2024 Engys Ltd.

Class
    Foam::discreteMixingPlanePolyPatch

Description
    This is a generalised frozen rotor approach.
    The complete 360 from both sides and all the possible sector-based
    combinations are constructed.
    The "classic" mixingPlane finds the average on different circumference
    levels/bands. In this approach the sectors are rotated in a discrete way
    depending on the sectors. If you have for example a 10-6 rotor stator
    sectors, the rotor will be rotated 10 times and the stator 6 times.
    All the possible combinations will be 60.
    Thus, it does not depend on the sector initial position as the the frozen
    rotor  approach (frozen rotor in literature).
    It has a small dependency.
    It is called “discrete” because discrete type of averaging is applied
    and not the band-based averaging as the classic mixingPlane.

Usage
    -sectorDefinition: periodicPatch
        Creates the averaging snapshots
        based on the periodic-side rotational transformations.
        Suggested for large number of sectors in both sides

        patch
        {
            type discreteMixingPlane;
            ...
            periodicPatch   cyclic01; // if there is
            sectorDefinition periodicPatch; //Default value
        }


    -sectorDefinition: userDefined
        Creates the averaging snapshots based on user input. Suggested for
        large sectors (or small nSector). Additional entries are needed to
        rotate the snapshots of each side and find the topological connections
        If the mesh has 3 sectors (fan) and the user wants to use x2
        snapshots then the nSectors must be 6 = 3(sectors)*2(nSnapshots)
        If the geometry is 360 then the nSectors represents the averaging
        snapshots.

        patch
        {
            type discreteMixingPlane;
            ...
            sectorDefinition userDefined;
            //- additional information needed for the transformation
            nSectors        16; //- 16 averaging snapshots
            sectorRotationAxis ( 0 0 1);
            sectorRotationCentre ( 0 0 0);
         }


SourceFiles
    discreteMixingPlanePolyPatch.C
    discreteMixingPlanePolyPatchSectors.C

\*---------------------------------------------------------------------------*/

#ifndef discreteMixingPlanePolyPatch_H
#define discreteMixingPlanePolyPatch_H

#include "AMIInterpolation/patches/cyclicAMI/cyclicAMIPolyPatch/cyclicAMIPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class OBJstream;

/*---------------------------------------------------------------------------*\
                    Class discreteMixingPlanePolyPatch Declaration
\*---------------------------------------------------------------------------*/

class discreteMixingPlanePolyPatch
:
    public cyclicAMIPolyPatch
{
public:

    // Public data types
        enum sectorDefinitions
        {
            periodicPatch,
            userDefined
        };

        static const NamedEnum<sectorDefinitions, 2> sectorDefinitionsTypes_;


private:

    // Private data

        //- Periodic patch name
        sectorDefinitions sectorDefinition_;

        //- Periodic patch name
        mutable word periodicPatchName_;

        //- Periodic patch ID
        mutable label periodicPatchID_;

        //- Ability to use different matching threshold for integer
        // sector checking in bad meshes and not the tolerance for the AMI
        const scalar sectorMatchTolerance_;

        //- Number of sectors in a rotationally periodic geometry
        //  Calculated automatically from periodics
        //  or user defines it for more averaging snapshots
        mutable label nSectors_;

        //- Number of subdivisions per sector
        //  Total snapshots will be nSectors_*nSubdivisions_
        mutable label nSubdivisions_;

        //- Transformations for each rotation
        mutable List<transformer> srcAMITransforms_;
        mutable List<transformer> tgtAMITransforms_;

        //- If sector definitions is provide by the user, the following two
        // entries are needed to compute the transformation
        mutable vector sectorRotationAxis_;
        mutable point sectorRotationCentre_;

        //- transformer of this and nbr patches
        mutable transformer thisTransformer_ = transformer::null;
        mutable transformer nbrTransformer_ = transformer::null;

    // Private Member Functions

        //- Check if 360 (periodicPatchName is null or none)
        bool is360() const;

        //- Synchronise the periodic transformations
        void syncTransforms() const;

        //- Compute nSectors
        void computeSector() const;

        //- Calculate and return transformer of this nbrghbour patch
        const transformer& nbrPatchTranformer() const;

        transformer getTransform(const coupledPolyPatch&) const;

        //- Debug: write obj files of patch (collected on master)
        void writeOBJ(const primitivePatch& p, OBJstream& str) const;

        //- Reset the AMI interpolator
        virtual void resetAMI
        (
            const AMIInterpolation::interpolationMethod& AMIMethod =
                AMIInterpolation::imFaceAreaWeight
        ) const;

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Number of sectrors
        virtual label nSectors() const;

        //- Periodic patch ID
        virtual label periodicPatchID() const;

        //- Neighbour periodic patch ID
        virtual label nbrPeriodicPatchID() const;

        //- Number of averaging snapshots
        virtual label thisSectors() const;

        //- Number of averaging snapshots of the nbrbour patch
        virtual label nbrSectors() const;



public:

    //- Runtime type information
    TypeName("discreteMixingPlane");


    // Constructors

        //- Construct from (base coupled patch) components
        discreteMixingPlanePolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from dictionary
        discreteMixingPlanePolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        discreteMixingPlanePolyPatch
        (
            const discreteMixingPlanePolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        discreteMixingPlanePolyPatch
        (
            const discreteMixingPlanePolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& nbrPatchName
        );

        //- Construct given the original patch and a map
        discreteMixingPlanePolyPatch
        (
            const discreteMixingPlanePolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );


        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>
            (
                new discreteMixingPlanePolyPatch(*this, bm)
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<directPolyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<directPolyPatch>
            (
                new discreteMixingPlanePolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    nbrPatchName_
                )
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<directPolyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        ) const
        {
            return autoPtr<directPolyPatch>
            (
                new discreteMixingPlanePolyPatch
                (
                    *this,
                    bm,
                    index,
                    mapAddressing,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~discreteMixingPlanePolyPatch(){};

    // Member Functions
        //- Calculate and return transformer of this patch
        const transformer& thisPatchTranformer() const;

        virtual bool implicitNonCovered() const
        {
            return false;
        }

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
