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
    (c) 2011-2018 OpenFOAM Foundation
    (c) 2010-2019 Engys Ltd.

Class
    Foam::mappedFvPatch

Description
    A region-coupled patch, which acts as an lduInterface for
    implicit coupling if used with a monolithicMatrix

SourceFiles
    mappedFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef mappedFvPatch_H
#define mappedFvPatch_H

#include "fvMesh/fvMesh.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyPatches/polyPatch/polyPatch.H"
#include "mappedPatches/mappedPolyPatch/mappedPolyPatch.H"
#include "mappedPatches/mappedLduInterface/mappedLduInterface.H"
#include "regionCoupled/patches/regionCoupledPolyPatch/regionCoupledPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class mappedFvPatch Declaration
\*---------------------------------------------------------------------------*/

class mappedFvPatch
:
    virtual public fvPatch,
    public mappedLduInterface,
    public lduInterface
{
    // Private data

        //- Reference to mappedPolyPatch
        const mappedPatchBase& mappedPolyPatch_;


protected:

    // Protected members

        //- Returns fvMesh
        const fvMesh& nbrFvMesh() const
        {
            return refCast<const fvMesh>(nbrMesh());
        }

        virtual void initMovePoints() override
        {
            // If this mesh moves but the other doesn't, the other
            // polyPatch does not clear its data as it needs to, e.g. to
            // recompute the AMI if it owns the AMI. Hence, force clear.
            const_cast<mappedPatchBase&>
            (
                dynamicCast<const mappedPatchBase>
                (
                    mappedPolyPatch_.samplePolyPatch()
                )
            ).clearOut();

            fvPatch::initMovePoints();
        }

public:

    //- Runtime type information
    TypeName(mappedPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        mappedFvPatch
        (
            const polyPatch& pp,
            const fvBoundaryMesh& bm
        )
        :
            fvPatch(pp, bm),
            mappedLduInterface(),
            mappedPolyPatch_
            (
                refCast<const mappedPatchBase>(pp)
            )
        {}


    //- Destructor
    virtual ~mappedFvPatch()
    {}


    // Member Functions

        // Access

            //- Whether it is able to act as an implicitly coupled interface
            virtual bool isInterface() const override
            {
                return
                (
                    (
                        this->mappedPolyPatch_.mode()
                        == mappedPolyPatch::NEARESTPATCHFACE
                     || this->mappedPolyPatch_.mode()
                        == mappedPolyPatch::NEARESTPATCHFACEAMI
                    )
                 && nbrMeshFound()
                 && isA<mappedFvPatch>(nbrFvMesh().boundary()[nbrPatchID()])
                );
            }

            //- Return neighbour
            virtual label nbrPatchID() const
            {
                return mappedPolyPatch_.samplePolyPatch().index();
            }

            //- Return mappedFvPatch neighbour patch
            virtual const mappedFvPatch& nbrPatch() const override
            {
                return refCast<const mappedFvPatch>
                (
                    nbrFvMesh().boundary()[nbrPatchID()]
                );
            }

            //- Does neighbour mesh exist?
            virtual bool nbrMeshFound() const override
            {
                return mappedPolyPatch_.sampleMeshFound();
            }

            //- Returns neighbour polyMesh
            virtual const polyMesh& nbrMesh() const override
            {
                return mappedPolyPatch_.sampleMesh();
            }

            //- Return fvPatch
            const fvPatch& patch() const
            {
                return *this;
            }

            //- Is it the owner of the AMI?
            virtual bool owner() const override
            {
                if (mappedPolyPatch_.sameRegion())
                {
                    return this->index() < nbrPatchID();
                }
                return
                (
                    this->boundaryMesh().mesh().name()
                  < this->nbrFvMesh().name()
                );
            }

            //- Return a reference to the AMI interpolator
            virtual const AMIInterpolation& AMI() const override
            {
                return mappedPolyPatch_.AMI();
            }

            //- Return reference to direct map
            virtual const distributionMap& map() const override
            {
                return mappedPolyPatch_.map();
            }

            //- Whether we are using AMI or direct map
            virtual bool usingAMI() const override
            {
                return
                (
                    mappedPolyPatch_.mode()
                 == mappedPolyPatch::NEARESTPATCHFACEAMI
                );
            }


            //- Return faceCell addressing
            virtual const labelUList& faceCells() const override
            {
                return fvPatch::faceCells();
            }

            //- Returns false - no interpolation across patch
            virtual bool coupled() const override
            {
                return false;
            }

            virtual const word& patchName() const override
            {
                return this->name();
            }

            virtual const word& regionName() const override
            {
                return this->boundaryMesh().mesh().name();
            }

            virtual const word& nbrPatchName() const override
            {
                return nbrPatch().name();
            }

            virtual const word& nbrRegionName() const override
            {
                return nbrPatch().boundaryMesh().mesh().name();
            }

            //- Return the interface type
            const word& regionCoupleType() const
            {
                return patch().type();
            }

            //- Name of the sample/neighbour region
            const word& sampleRegion() const
            {
                return mappedPolyPatch_.sampleRegion();
            }

            //- Map or use the AMI to perform interpolation from neighbour
            // defaultValues are used for AMI low-weight correction
            template <class Type>
            tmp<Field<Type>> interpolate
            (
                const Field<Type>& fld,
                const UList<Type>& defaultValues
            ) const;

            //- Map or use the AMI to perform interpolation from neighbour
            // defaultValues are used for AMI low-weight correction
            template <class Type>
            tmp<Field<Type>> interpolate
            (
                tmp<Field<Type>> tfld,
                const UList<Type>& defaultValues
            ) const;

    // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const override;

            //- Initialise neighbour field transfer
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const override
            {}


            //- Return neighbour field. Called when the solve is being done
            //  segregated
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const override;

            //- Return neighbour field. Called when the solve is being
            //  done region-coupled
            virtual tmp<labelField> regionCoupledInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvMesh/fvPatches/derived/mapped/mappedFvPatchTemplates.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
