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
    (c) 2015-2016 OpenFOAM Foundation
    (c) 2022 Engys Ltd.

Class
    Foam::distributedWeightedFvPatchFieldMapper

Description
    FieldMapper with weighted mapping from (optionally remote) quantities.

\*---------------------------------------------------------------------------*/

#ifndef distributedWeightedFvPatchFieldMapper_H
#define distributedWeightedFvPatchFieldMapper_H

#include "fields/fvPatchFields/fvPatchField/fieldMappers/distributedFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class distributedWeightedFvPatchFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class distributedWeightedFvPatchFieldMapper
:
    public distributedFvPatchFieldMapper
{
    // Private Data

        const label singlePatchProc_;

        const labelListList& addressing_;

        const scalarListList& weights_;


    // Private Member Functions

        label size() const
        {
            if (distributed())
            {
                return distMapPtr_->constructSize();
            }
            else
            {
                return addressing().size();
            }
        }

        template<class Type>
        void map
        (
            Field<Type>& f,
            const Field<Type>& mapF,
            const bool applyFlip
        ) const;

        template<class Type>
        tmp<Field<Type>> map
        (
            const Field<Type>& f,
            const bool applyFlip
        ) const;


public:

    // Constructors

        //- Construct given addressing
        distributedWeightedFvPatchFieldMapper
        (
            const label singlePatchProc,
            const distributionMapBase* distMapPtr,
            const labelListList& addressing,
            const scalarListList& weights
        )
        :
            distributedFvPatchFieldMapper(distMapPtr),
            singlePatchProc_(singlePatchProc),
            addressing_(addressing),
            weights_(weights)
        {
            forAll(addressing_, i)
            {
                if (addressing_[i].size() == 0)
                {
                    hasUnmapped_ = true;
                }
            }

            if ((singlePatchProc_ == -1) != (distMapPtr_ != nullptr))
            {
                FatalErrorInFunction
                    << "Supply a distributionMapBase if and only if "
                    << "singlePatchProc is -1"
                    << " singlePatchProc_:" << singlePatchProc_
                    << " distMapPtr_:" << (distMapPtr_ != nullptr)
                    << exit(FatalError);
            }
        }


    //- Destructor
    virtual ~distributedWeightedFvPatchFieldMapper()
    {}


    // Member Functions

        //- Is the mapping direct
        virtual bool direct() const override
        {
            return false;
        }

        //- No indirect mapping
        virtual bool indirect() const override
        {
            return false;
        }

        //- Is the mapping distributed
        virtual bool distributed() const override
        {
            return singlePatchProc_ == -1;
        }

        //- Are there unmapped values? I.e. do all size() elements get value?
        virtual bool hasUnmapped() const override
        {
            return hasUnmapped_;
        }

        //- Return the interpolated map addressing
        virtual const labelListList& addressing() const override
        {
            return addressing_;
        }

        //- Return interpolation weights
        virtual const scalarListList& weights() const override
        {
            return weights_;
        }


    // Member Operators

        //- Map a field
        FOR_ALL_FIELD_TYPES(DECLARE_DISTRIBUTED_MAPPER_OPERATOR, , override)

        //- Map a label field
        DECLARE_DISTRIBUTED_MAPPER_OPERATOR(label, , override)

        //- Specialisations for VectorN types
        forAllVectorNTypes(DECLARE_DISTRIBUTED_MAPPER_OPERATOR, override)

        forAllTensorNTypes(DECLARE_DISTRIBUTED_MAPPER_OPERATOR, override)

        forAllDiagTensorNTypes(DECLARE_DISTRIBUTED_MAPPER_OPERATOR, override)

        forAllSphericalTensorNTypes
        (
            DECLARE_DISTRIBUTED_MAPPER_OPERATOR,
            override
        )
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "meshToMesh/distributedWeightedFvPatchFieldMapperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
