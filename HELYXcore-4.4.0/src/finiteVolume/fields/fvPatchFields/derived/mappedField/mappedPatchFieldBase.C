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
    (c) 2016-2024 Engys Ltd
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/mappedField/mappedPatchFieldBase.H"
#include "mappedPatches/mappedPolyPatch/mappedPatchBase.H"
#include "interpolation/interpolation/interpolationCell/interpolationCell.H"
#include "fields/Fields/symmTransformField/symmTransformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
mappedPatchFieldBase<Type>::mappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField,
    const word& fieldName,
    const bool setAverage,
    const Type average,
    const word& interpolationScheme
)
:
    mapper_(mapper),
    patchField_(patchField),
    fieldName_(fieldName),
    setAverage_(setAverage),
    average_(average),
    flowRate_(Zero),
    isMassFlow_(true),
    rhoName_("rho"),
    interpolationScheme_(interpolationScheme)
{}


template<class Type>
mappedPatchFieldBase<Type>::mappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField,
    const dictionary& dict
)
:
    mapper_(mapper),
    patchField_(patchField),
    fieldName_
    (
        dict.template lookupOrDefault<word>
        (
            "field",
            patchField_.internalField().name()
        )
    ),
    setAverage_(dict.lookup<bool>("setAverage")),
    average_
    (
        ( setAverage_ && dict.found("average") )
        ? pTraits<Type>(dict.lookup("average"))
        : dict.lookupOrDefault<Type>("average", Zero)
    ),
    flowRate_
    (
        ( setAverage_ && dict.found("flowRate") )
        ? dict.lookup<scalar>("flowRate")
        : 0.0
    ),
    isMassFlow_(dict.lookupOrDefault<Switch>("isMassFlowRate", true)),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    interpolationScheme_(interpolationCell<Type>::typeName)
{
    if (mapper_.mode() == mappedPatchBase::NEARESTCELL)
    {
        interpolationScheme_ = dict.lookup<word>("interpolationScheme");
    }
}


template<class Type>
mappedPatchFieldBase<Type>::mappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField
)
:
    mapper_(mapper),
    patchField_(patchField),
    fieldName_(patchField_.internalField().name()),
    setAverage_(false),
    average_(Zero),
    flowRate_(Zero),
    isMassFlow_(true),
    rhoName_("rho"),
    interpolationScheme_(interpolationCell<Type>::typeName)
{}


template<class Type>
mappedPatchFieldBase<Type>::mappedPatchFieldBase
(
    const mappedPatchFieldBase<Type>& mapper
)
:
    mapper_(mapper.mapper_),
    patchField_(mapper.patchField_),
    fieldName_(mapper.fieldName_),
    setAverage_(mapper.setAverage_),
    average_(mapper.average_),
    flowRate_(mapper.flowRate_),
    isMassFlow_(mapper.isMassFlow_),
    rhoName_(mapper.rhoName_),
    interpolationScheme_(mapper.interpolationScheme_)
{}


template<class Type>
mappedPatchFieldBase<Type>::mappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField,
    const mappedPatchFieldBase<Type>& base
)
:
    mapper_(mapper),
    patchField_(patchField),
    fieldName_(base.fieldName_),
    setAverage_(base.setAverage_),
    average_(base.average_),
    flowRate_(base.flowRate_),
    isMassFlow_(base.isMassFlow_),
    rhoName_(base.rhoName_),
    interpolationScheme_(base.interpolationScheme_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const VolField<Type>&
mappedPatchFieldBase<Type>::sampleField() const
{
    typedef VolField<Type> fieldType;

    const fvMesh& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());

    if (mapper_.sameRegion())
    {
        if (fieldName_ == patchField_.internalField().name())
        {
            // Optimisation: bypass field lookup
            return
                dynamic_cast<const fieldType&>
                (
                    patchField_.internalField()
                );
        }
        else
        {
            const objectRegistry& thisReg = patchField_.db();
            return thisReg.template lookupObject<fieldType>(fieldName_);
        }
    }
    else
    {
        return nbrMesh.template lookupObject<fieldType>(fieldName_);
    }
}


template<class Type>
Foam::scalar mappedPatchFieldBase<Type>::getUfromFlowRate() const
{
    FatalErrorInFunction
        << "Function is not implemented in the base class"
        << exit(FatalError);

    return Zero;
}


template<class Type>
tmp<Field<Type>> mappedPatchFieldBase<Type>::mappedField() const
{
    typedef VolField<Type> fieldType;

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    const fvMesh& thisMesh = patchField_.patch().boundaryMesh().mesh();
    const fvMesh& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());

    // Result of obtaining remote values
    tmp<Field<Type>> tnewValues(new Field<Type>(0));
    Field<Type>& newValues = tnewValues.ref();

    switch (mapper_.mode())
    {
        case mappedPatchBase::NEARESTCELL:
        {
            const distributionMap& distMap = mapper_.map();

            if (interpolationScheme_ != interpolationCell<Type>::typeName)
            {
                // Send back sample points to the processor that holds the cell
                vectorField samples(mapper_.samplePoints());
                distMap.reverseDistribute
                (
                    (
                        mapper_.sameRegion()
                      ? thisMesh.nCells()
                      : nbrMesh.nCells()
                    ),
                    point::max,
                    samples
                );

                autoPtr<interpolation<Type>> interpolator
                (
                    interpolation<Type>::New
                    (
                        interpolationScheme_,
                        sampleField()
                    )
                );
                const interpolation<Type>& interp = interpolator();

                newValues.setSize(samples.size(), pTraits<Type>::max);
                forAll(samples, celli)
                {
                    if (samples[celli] != point::max)
                    {
                        newValues[celli] = interp.interpolate
                        (
                            samples[celli],
                            celli
                        );
                    }
                }
            }
            else
            {
                newValues = sampleField();
            }

            distMap.distribute(newValues);

            break;
        }
        case mappedPatchBase::NEARESTPATCHFACE:
        case mappedPatchBase::NEARESTPATCHFACEAMI:
        {
            const label nbrPatchID =
                nbrMesh.boundaryMesh().findPatchID(mapper_.samplePatch());

            if (nbrPatchID < 0)
            {
                FatalErrorInFunction
                 << "Unable to find sample patch " << mapper_.samplePatch()
                 << " in region " << mapper_.sampleRegion()
                 << " for patch " << patchField_.patch().name() << nl
                 << abort(FatalError);
            }

            const fieldType& nbrField = sampleField();

            newValues = nbrField.boundaryField()[nbrPatchID];
            mapper_.distribute
            (
                newValues,
                dynamic_cast<const UList<Type>&>
                (
                    symmetryValue
                    (
                        patchField_.patchInternalField(),
                        patchField_.patch().Sf()
                    )()
                )
            );

            break;
        }
        case mappedPatchBase::NEARESTFACE:
        {
            Field<Type> allValues(nbrMesh.nFaces(), Zero);

            const fieldType& nbrField = sampleField();

            forAll(nbrField.boundaryField(), patchi)
            {
                const fvPatchField<Type>& pf =
                    nbrField.boundaryField()[patchi];
                label faceStart = pf.patch().start();

                forAll(pf, facei)
                {
                    allValues[faceStart++] = pf[facei];
                }
            }

            mapper_.distribute(allValues);
            newValues.transfer(allValues);

            break;
        }
        default:
        {
            FatalErrorInFunction
             << "Unknown sampling mode: " << mapper_.mode()
             << nl << abort(FatalError);
        }
    }

    scaleWithAverage(newValues);

    // Restore tag
    UPstream::msgType() = oldTag;

    return tnewValues;
}


template<class Type>
void mappedPatchFieldBase<Type>::scaleWithAverage(Field<Type>& fieldValue) const
{
    if (setAverage_)
    {
        Type averagePsi =
            gSum(patchField_.patch().magSf()*fieldValue)
            /gSum(patchField_.patch().magSf());

        if (mag(averagePsi)/mag(average_) > 0.5)
        {
            fieldValue *= mag(average_)/mag(averagePsi);
        }
        else
        {
            fieldValue += (average_ - averagePsi);
        }
    }
}


template<class Type>
void mappedPatchFieldBase<Type>::write(Ostream& os) const
{
    os.writeEntry("field", fieldName_);
    os.writeEntry("setAverage", setAverage_);
    os.writeEntry("average", average_);
    os.writeEntry("flowRate", flowRate_);
    os.writeEntry("isMassFlowRate", isMassFlow_);
    os.writeEntry("interpolationScheme", interpolationScheme_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
